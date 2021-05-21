import unittest
import os
import glob
import csv

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
                    truth_data = list( csv.reader(f, delimiter=' ') )
                    
                actual_file = os.path.join(actual_dir, ifile)
                with open(actual_file) as f:
                    actual_data = list( csv.reader(f, delimiter=' ') )

                # Print progress
                print(f'... now testing {truth_file}')

                # Loop over all file contents and compare one value at a time
                for iline in range(len(truth_data)):
                    for ikey in range(len(truth_data[iline])):
                        if isfloat(truth_data[iline][ikey]):
                            truth_float = float( truth_data[iline][ikey] )
                            actual_float = float( actual_data[iline][ikey] )
                            self.assertAlmostEqual(truth_float, actual_float, 4)
                        else:
                            self.assertEqual(truth_data[iline][ikey], actual_data[iline][ikey])

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
