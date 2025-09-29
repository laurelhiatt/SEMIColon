#mamba install SigProfilerMatrixGenerator

from SigProfilerMatrixGenerator import install as genInstall
genInstall.install('GRCh38', rsync=False, bash=True)


from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen
matrices = matGen.SigProfilerMatrixGeneratorFunc("test", "GRCh38", "/Users/ebergstr/Desktop/test",plot=True, exome=False, bed_file=None, chrom_based=False, tsb_stat=True, seqInfo=True, cushion=100)