To use the Tesseract.f90 code to perform the word ordering task:

1. Use the Makefile to "make tess".

2. Set the control file tesseract.txt to the desired parameters and pathnames (example control file provided).

3. Run HHM on a corpus (see https://github.com/MatthewAKelly/BEAGLE-HHM for MatLab Code or https://github.com/moojan/Python-BEAGLE-HHM for Python code).

4. test_file2017.txt is the set of test sentences used to evaluate the simple exemplar model in the word ordering task.

5. Give tess the memory vectors from the desired level of HHM to represent words and test_file2017.txt in order to evaluate the ability of the vectors to order words into sentences using a simple exemplar model.

To use beagle_pos_5.py to classify words by their part of speech:

1. Make sure that the part of speech category file, words.txt, is in the same folder.

2. Set the correct path names for the csv word embedding files in the code.

3. run beagle_pos_5.py in Python.