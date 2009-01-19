import unittest
import os
import order_boundary as ob


class Test_order_boundary(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_Simple(self):
        import os
        import csv

        # filenames
        Test_input_file_path = 'test.in.csv'
        Test_output_file_path = 'test.out.csv'

        # input data
        Data = [('longitude','latitude','index'),
                ( 1.0,  1.0, 'alpha', 'extra'),
                ( 2.3,  2.0, 'bravo'),
                ( 3.9,  3.0, 'charlie'),
                ( 9.0,  9.9, 'delta'),
                (10.0, 10.4, 'echo'),
                (11.0, 11.0, 'foxtrot'),
                (15.9, 16.0, 'golf'),
                (17.0, 17.1, 'hotel'),
                (17.9, 18.0, 'india'),
                (12.2, 12.0, 'juliet'),
                ( 4.7,  4.0, 'kilo'),
                ( 5.2,  5.0, 'lima'),
                ( 6.0,  6.0, 'mike'),
                ( 7.0,  7.3, 'november', 'extra', 'extra', 'extra'),
                ( 8.0,  8.7, 'oscar'),
                (13.6, 13.0, 'papa'),
                (14.9, 14.0, 'quebec'),
                (15.8, 15.0, 'romeo'),
                (16.0, 16.2, 'sierra'),
                (17.1, 17.1, 'tango'),
                (18.0, 18.7, 'uniform'),
                (19.0, 19.9, 'victor'),
                (20.0, 20.0, 'whisky')]

        # expected output data
        Expected = [('longitude','latitude','index'),
                    ( 1.0,  1.0, 'alpha', 'extra'),
                    ( 2.3,  2.0, 'bravo'),
                    ( 3.9,  3.0, 'charlie'),
                    ( 4.7,  4.0, 'kilo'),
                    ( 5.2,  5.0, 'lima'),
                    ( 6.0,  6.0, 'mike'),
                    ( 7.0,  7.3, 'november', 'extra', 'extra', 'extra'),
                    ( 8.0,  8.7, 'oscar'),
                    ( 9.0,  9.9, 'delta'),
                    (10.0, 10.4, 'echo'),
                    (11.0, 11.0, 'foxtrot'),
                    (12.2, 12.0, 'juliet'),
                    (13.6, 13.0, 'papa'),
                    (14.9, 14.0, 'quebec'),
                    (15.8, 15.0, 'romeo'),
                    (15.9, 16.0, 'golf'),
                    (16.0, 16.2, 'sierra'),
                    (17.0, 17.1, 'hotel'),
                    (17.1, 17.1, 'tango'),
                    (17.9, 18.0, 'india'),
                    (18.0, 18.7, 'uniform'),
                    (19.0, 19.9, 'victor'),
                    (20.0, 20.0, 'whisky')]

        # put test data into a file
        fd = open(Test_input_file_path, 'wb')
        w = csv.writer(fd)
        for d in Data:
            w.writerow(d)
        fd.close()

        # call routine, put sorted points into output file
        ob.order_boundary(Test_input_file_path, Test_output_file_path)

        # get sorted data into memory
        fd = open(Test_output_file_path, 'r')
        data_list = []
        for data in csv.reader(fd):
            try:
                data[0] = float(data[0])
            except:
                pass
            try:
                data[1] = float(data[1])
            except:
                pass
            data_list.append(tuple(data))
        fd.close()

        # check same as Expected
        self.failUnless(data_list == Expected)

        # clean up
        try:
            os.remove(Test_input_file_path)
        except:
            pass
        try:
            os.remove(Test_output_file_path)
        except:
            pass

    def test_OK(self):
        infile = 'test.csv'
        outfile = 'test.out.csv'

        ob.order_boundary(infile, outfile)

        os.remove(outfile)


#-------------------------------------------------------------
if __name__ == "__main__":
    suite = unittest.makeSuite(Test_order_boundary, 'test')
    runner = unittest.TextTestRunner()
    runner.run(suite)

