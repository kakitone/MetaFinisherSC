import unittest


class test(unittest.TestCase):
    def setUp(self):
        print "Init"
        
    def testone(self):
        print "1: "


def main():
    unittest.main()
    
if __name__ == '__main__':
    main()
