import os, glob
import subprocess
import argparse
from stefanHMM_concurrent import *
from hmmSearcher import * 
import os
import subprocess

# evol = 'M4'
# tmpdir = "/Users/gillianchu/warnow/bin/gitrepos/stefan_code_latest/alignData/tmpfiles/"
# datapath = ""/Users/gillianchu/warnow/bin/gitrepos/stefan_code_latest/alignData/UnalignFragTree/high_frag/1000" + evol + "/""
# bundle_packagedir = "/Users/gillianchu/.sepp/bundled-v4.3.17/"
# packagedir = ""
# outdir =  "/Users/gillianchu/warnow/bin/gitrepos/stefan_code_latest/alignData/UPPoutput/"

def parse_strats(stratsfile): 
    stratlist = []
    with open(stratsfile, "r") as reader:
        lines = reader.readlines()
        for line in lines:
            stratlist.append(line.strip())
    return stratlist

def build_upp_config(evol, root_tmpdir, indata, bundle_packagedir, packagedir, numreps, decomp_size, outdir):
    if not os.path.exists(root_tmpdir + evol):
        subprocess.call(['mkdir', root_tmpdir + evol])

    for r in range(numreps):    
        rep = "R"+str(r)
        unaligned_frag_dir = indata + rep + "/unaligned_frag.txt"
        unaligned_frag = glob.glob(unaligned_frag_dir)[0] 
        assert len(glob.glob(unaligned_frag_dir)) == 1
        print("Running...", unaligned_frag, evol + rep + "_decomp" + str(decomp_size))
        
        currdir = indata + rep + "/"

        if not os.path.exists(root_tmpdir + evol + "/" + rep):
            subprocess.call(['mkdir', root_tmpdir + evol + "/" + rep])

        tmpdir = root_tmpdir + evol + "/" + rep + "/" + str(decomp_size)

        alignfile = currdir + "pasta_backbone_align.txt"
        treefile = currdir + "pasta_backbone.tre"

        if not os.path.exists(tmpdir):
            subprocess.call(['mkdir', tmpdir])

        configfile = root_tmpdir + evol + "/" + rep + "/" + "decomp" + str(decomp_size) + ".config"
        with open(configfile, "w+") as f:
            f.write("[pplacer]\n")
            f.write("path=%spplacer\n" % bundle_packagedir)
            f.write("[hmmalign]\n")
            f.write("path=%shmmalign\n" % packagedir)
            f.write("[hmmsearch]\npath=%shmmsearch\npiped=False\nelim=10000\nfilters=True\n" % packagedir)
            f.write("[hmmbuild]\npath=%shmmbuild\n" % packagedir)
            f.write("[jsonmerger]\npath=%sseppJsonMerger.jar\n" % bundle_packagedir)
            f.write("[exhaustive]\nstrategy = centroid\nminsubsetsize = %s\nplacementminsubsetsizefacotr = 4\nplacer = pplacer\nweight_placement_by_alignment = True\n" % str(decomp_size))
            f.write('[commandline]\n')
            f.write("alignment=%s\ntree=%s\nsequence_file=%s\n"%(alignfile,treefile,unaligned_frag))
            f.write("alignmentSize=%s\n" % str(decomp_size))
            f.write("molecule=dna\ntempdir=%s\n" % tmpdir)
            f.write("cpu=1")

        subprocess.call(['python', 'run_upp.py', "-c", configfile])

        outputfile = outdir + evol + rep + "decomp" + str(decomp_size)
        subprocess.call(['mv', 'output_alignment.fasta', outputfile + "_output_alignment" + ".fasta"])
        subprocess.call(['mv', 'output_alignment_masked.fasta', outputfile + "_output_alignment_masked.fasta"])
        subprocess.call(['mv', 'output_insertion_columns.txt', outputfile + "_output_insertion_columns.txt"])

def makedirs(dirpath): 
    # clear the folders
    for firstlvl in ["trueAlignment","temporaryFileSave","subsetTrueAln", "Searcher","queryToHmm","newHMM","ML","initialHMM","hmmSeqAlign","hmmScores","hmmQueryList","fullPrediction"]:
        
        firstdir = '%s/%s' % (dirpath, firstlvl)
        if os.path.isdir(firstdir):
            subprocess.call(['rm', firstdir])

        subprocess.call(['mkdir', firstdir])
        if firstlvl == 'fullPrediction':
            subprocess.call(['mkdir', firstdir + '/sequences/'])
        elif firstlvl == 'hmmQueryList':
            for z in ['inputQuery', 'merged', 'predictedQuery']:
                subprocess.call(['mkdir', firstdir + '/' + z])
        elif firstlvl == 'hmmScores':
            for z in ['newRawOld', 'processedOld', 'rawOld', 'scoresFull', 'temporaryStorage']:
                subprocess.call(['mkdir', firstdir + '/' + z])
        elif firstlvl == 'ML':
            for z in ['hmm', 'models', 'queryHMM', 'scores', 'sequences', 'temporaryStorage']:
                subprocess.call(['mkdir', firstdir + '/' + z])
            subprocess.call(['mkdir', firstdir + '/sequences/np'])
        elif firstlvl == 'newHMM':
            for z in ['columnSets', 'hmm', 'newHMMseq']:
                subprocess.call(['mkdir', firstdir + '/' + z])
        elif firstlvl == 'queryToHmm':
            for z in ['original', 'withUPP']:
                subprocess.call(['mkdir', firstdir + '/' + z])
        elif firstlvl == 'Searcher':
            subprocess.call(['mkdir', '%s/scoreFiles' % firstdir])
        elif firstlvl == 'trueAlignment':
            for z in ['original', 'subset']:
                subprocess.call(['mkdir', firstdir + '/' + z])
       
def run_upp_strats(evol, dirpath, numreps, strats, doResort=False):

    for i in range(numreps):  
        rep = "R%s" % str(i)
        decomp = str(d)

        print("[START HERE] evol:%s, rep:%s, decomp:%s" % (evol, rep, decomp), flush=True)
        dirname = "./%s/tmpfiles/%s/%s/%s/" % (dirpath, evol, rep, decomp)
        outputdirname = os.listdir(dirname)[0]
        hmmSeqFile = '%s/%s/root/P_0/' % (dirname, outputdirname)
        fragmentfile = os.listdir(dirname + outputdirname + "/fragment_chunks/")[0]
        queryName = "%s/%s/fragment_chunks/%s" % (dirname, outputdirname, fragmentfile)

        trueAlignment = "./%s/UnalignFragTree/high_frag/1000%s/%s/true_align_fragged.txt" % (dirpath, evol, rep)
        predictionName = './%s/UPPoutput/%s%sdecomp%s_output_alignment.fasta' % (dirpath, evol, rep, decomp)

        setAllFileNames(hmmSeqFile, queryName, trueAlignment, predictionName)
        saveInitialSteps()
        ## saveDecomposition() new version contains this

        hierchySearch()

        for strat in strats: 
            if strat != 'stefan_trueUPP':
                print("[processing %s]" % strat)
                print("[running scoresToHMMSeq]")
                scoresToHMMSeq(strat)
                print("[running buildAlignMerge, doResort is %s]" % doResort)
                buildAlignMerge(strat, doResort=doResort)
            print("[running scoreAlignment]")
            scoreAlignment(strat)

def main(): 
    parser = argparse.ArgumentParser()

    ## take in upp tree
    ## take in upp backbone alignment

    parser.add_argument('evol', type=str, help="evolutionary model i.e. M1 or M4", required=True)
    parser.add_argument('tmpdir', type=str, help="absolute path for UPP tmp files", required=True)
    parser.add_argument('datapath', type=str, help="absolute path for data", required=True)
    parser.add_argument('hmmer_packagedir', type=str, help="hmmer package dir", required=True)
    parser.add_argument('bundle_packagedir', type=str, help="bundled package dir", required=True)
    parser.add_argument('numreps', type=int, help="num reps", required=True)
    parser.add_argument('method', type=str, help="see README for more details on each method", required=True)
    parser.add_argument('decomp', type=int, help="decomp minimumsubset size")
    parser.add_argument('doResort', type=str, help="True/False, default is False")
    parser.add_argument('strats', type=str, help='file with strategies to run one per line, see README for more details')
    args = parser.parse_args()

    evol = args.evol
    root_tmpdir = args.tmpdir 
    indata = args.datapath
    bundle_packagedir = args.bundle_packagedir
    packagedir = args.hmmer_packagedir
    numreps = args.numreps
    decomp = args.decomp
    doResort = args.doResort
    strats = parse_strats(args.strats)

    dirpath = 'ensemble_data'
    outdir = '%s/UPPOutput/' % dirpath 
    build_upp_config(evol, root_tmpdir, indata, bundle_packagedir, packagedir, numreps, decomp, dirpath)
    makedirs(dirpath)
    run_upp_strats(evol, dirpath, numreps, strats, doResort)

if __name__ == '__main__':
    main()
