import unittest
import os
import glob
import math

def isfloat(instr):
    try:
        _ = float(instr)
        return True
    except:
        return False

    
class TestCertRegression(unittest.TestCase):

    def testAll(self):
        certs = ['Cylinder','DeepCwind','HywindSpar','Moonpool']
        start_dir = os.path.dirname( os.path.realpath(__file__) )
        
        for icert in certs:
            print(f'Running regression tests for {icert} example')
            # Set truth and new sample directories
            truth_dir = os.path.join(start_dir, icert, 'Output_Benchmark')
            actual_dir = os.path.join(start_dir, icert, 'Output')

            # Get all files in the truth directory
            os.chdir(truth_dir)
            all_files = glob.glob(os.path.join('**', '*.*'), recursive=True)
            all_files.sort()
            
            # Loop over all files
            for ifile in all_files:

                # Load file contents into lists
                truth_file = os.path.join(truth_dir, ifile)
                with open(truth_file) as f:
                    truth_data = list( f.read().splitlines() )
                    
                actual_file = os.path.join(actual_dir, ifile)
                with open(actual_file) as f:
                    actual_data = list( f.read().splitlines() )

                # Print progress
                print(f'... now testing {truth_file}')

                # Loop over all file contents and compare one value at a time
                for iline in range(len(truth_data)):
                    truth_tok = truth_data[iline].split()
                    actual_tok = actual_data[iline].split()
                    for ikey in range(len(truth_tok)):
                        if isfloat(truth_tok[ikey]) and isfloat(actual_tok[ikey]):
                            truth_float = float( truth_tok[ikey] )
                            actual_float = float( actual_tok[ikey] )
                            answer = math.isclose(truth_float, actual_float, rel_tol=1e-1, abs_tol=0.0)
                            self.assertTrue(answer, f'TRUTH {truth_float} vs NEW {actual_float}')
                        elif ( (isfloat(truth_tok[ikey]) and not isfloat(actual_tok[ikey])) or
                               (not isfloat(truth_tok[ikey]) and isfloat(actual_tok[ikey])) ):
                               breakpoint()
                        else:
                            self.assertEqual(truth_tok[ikey], actual_tok[ikey])

            # Change back to starting point
            os.chdir(start_dir)

            

def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestCertRegression))
    return suite


if __name__ == "__main__":
    result = unittest.TextTestRunner().run(suite())

    if result.wasSuccessful():
        exit(0)
    else:
        exit(1)
