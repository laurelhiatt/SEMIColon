from SigProfilerExtractor import sigpro as sig
import sys

input_type = "matrix"
output = "ID"
input_data = "/uufs/chpc.utah.edu/common/HIPAA/u1264408/u1264408/Git/SEMIColon/data/output/CellCut/overlap_variants/indels/output/ID/indels.ID83.all"

def run():
    sig.sigProfilerExtractor(input_type,
                     output,
                     input_data,
                     reference_genome="GRCh38",
                     opportunity_genome = "GRCh38",
                     context_type = "ID",
                     exome = False,
                     minimum_signatures=1,
                     maximum_signatures=15,
                     nmf_replicates=100,
                     resample = True,
                     batch_size=1,
                     cpu=-1,
                     gpu=False,
                     nmf_init="random",
                     precision= "single",
                     matrix_normalization= "gmm",
                     seeds= "random",
                     min_nmf_iterations= 10000,
                     max_nmf_iterations=1000000,
                     nmf_test_conv= 10000,
                     nmf_tolerance= 1e-15,
                     get_all_signature_matrices= False)

if __name__ == '__main__':
    # sample = sys.argv[1]
    # input = './' + sample + '_in'
    # output = './' + sample +        '_out'
    run()