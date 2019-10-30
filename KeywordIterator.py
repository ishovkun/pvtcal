#!/usr/bin/env python

from WordIterator import WordIterator
from copy import copy

class SkipCommentIterator(WordIterator):
    """
    Same as the WordIterator class but skips
    line comments
    """
    comment_symbol = "--"

    def __next__(self):
        val = WordIterator.__next__(self)

        if (val[0:len(self.comment_symbol)] != self.comment_symbol):
            return val
        else:                   # skip line
            self.words = next(self.it_line).split()
            while (len(self.words) == 0):
                self.words = next(self.it_line).split()

            self.it_word = iter(self.words)
            return self.__next__()


class KeywordIterator:
    """
    Iterate over Eclipse File Keywords
    """
    word_it = None              # WordIterator
    terminal_char = "/"
    double_terminated_words = set()

    def __init__(self, fname):
        self.word_it = SkipCommentIterator(fname)
        self.initDoubleTerminatedSet_()

    def __iter__(self):
        return self

    def __next__(self):
        word = next(self.word_it)
        if (word in self.double_terminated_words):
            return word, self.readDoubleTerminatedKwd_()
        else:
            return word, self.readSingleTerminatedKwd_()
            
    def checkWordTerminates_(self, word):
        split = word.split(self.terminal_char)
        if (len(split) > 1 and len(split[-1]) == 0):
            return True
        elif (len(split) > 1 and len(split[-1]) > 0):
            raise ValueError("there should be a whitespace after / ")
        else: return False

    def readSingleTerminatedKwd_(self):
        words = []
        word = next(self.word_it)
        while ( not self.checkWordTerminates_(word) ):
            words.append(word)
            word = next(self.word_it)

        # do not return the "/" char
        self.handleTerminating(word, words)
        return words

    def handleTerminating(self, word, appendTo):
        assert( self.checkWordTerminates_(word) )
        spl = word.split(self.terminal_char)
        if (len(spl[0]) > 0):
            appendTo.append(spl[0])

    def readDoubleTerminatedKwd_(self):
        parts = []
        words = []
        interrupt = False
        while (not interrupt):
            word = next(self.word_it)
            while ( not self.checkWordTerminates_(word) ):
                words.append(word)
                word = next(self.word_it)
            self.handleTerminating(word, words)
            if (len(words) == 0):
                interrupt = True
            else:
                parts.append(copy(words))
                words = []

        return parts

    def initDoubleTerminatedSet_(self):
        self.double_terminated_words = set([
            "NONLINEAR", "WCONINJE", "WCONPROD", "WELLSTRE", # wells
            "SCOND", "PVTO"                                  # fluid
        ])
