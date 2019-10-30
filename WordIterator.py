#!/usr/bin/env python


class WordIterator:
    """
    This class implements the missing python iterator
    by word in a file
    """
    file = None       # current file object
    words = None      # list of current line's words splitted by whitespace
    it_line = None    # current line iterator (iterate over file)
    it_word = None    # cuttent word (iterate over words list)

    def __init__(self, file_name):
        self.file = open(file_name, "r")
        self.it_line = iter(self.file)
        self.words = next(self.it_line).split()
        self.it_word = iter(self.words)

    def __del__(self):
        self.file.close()
        it_line = None
        it_word = None

    def __iter__(self):
        return self

    def __next__(self):
        if (self.it_word is None):
            raise ValueError("The file has been closed")

        # read from current line
        try:
            return next(self.it_word)

        # take next line
        except StopIteration:
            self.words = next(self.it_line).split()
            while (len(self.words) == 0):
                self.words = next(self.it_line).split()

            self.it_word = iter(self.words)
            return next(self.it_word)
