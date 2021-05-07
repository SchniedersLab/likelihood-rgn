#!/usr/bin/python

# imports
import sys
import os
import re
import tensorflow as tf

# Constants
NUM_DIMENSIONS = 3

# Accesory functions for dealing with TF Example and SequenceExample
_example = tf.train.Example
_sequence_example = tf.train.SequenceExample
_feature = tf.train.Feature
_features = lambda d: tf.train.Features(feature=d)
_feature_list = lambda l: tf.train.FeatureList(feature=l)
_feature_lists = lambda d: tf.train.FeatureLists(feature_list=d)
_bytes_feature = lambda v: _feature(bytes_list=tf.train.BytesList(value=v))
_int64_feature = lambda v: _feature(int64_list=tf.train.Int64List(value=v))
_float_feature = lambda v: _feature(float_list=tf.train.FloatList(value=v))

# Functions for conversion from Text-based ProteinNet files to TFRecords
_aa_dict = {'A': '0', 'C': '1', 'D': '2', 'E': '3', 'F': '4', 'G': '5', 'H': '6', 'I': '7', 'K': '8', 'L': '9', 'M': '10', 'N': '11', 'P': '12', 'Q': '13', 'R': '14', 'S': '15', 'T': '16', 'V': '17', 'W': '18', 'Y': '19'}
_dssp_dict = {'L': '0', 'H': '1', 'B': '2', 'E': '3', 'G': '4', 'I': '5', 'T': '6', 'S': '7'}
_mask_dict = {'-': '0', '+': '1'}

class switch(object):
    """Switch statement for Python, based on recipe from Python Cookbook."""

    def __init__(self, value):
        self.value = value
        self.fall = False

    def __iter__(self):
        """Return the match method once, then stop"""
        yield self.match
        raise StopIteration
    
    def match(self, *args):
        """Indicate whether or not to enter a case suite"""
        if self.fall or not args:
            return True
        elif self.value in args: # changed for v1.5
            self.fall = True
            return True
        else:
            return False

def letter_to_num(string, dict_):
    """ Convert string of letters to list of ints """
    patt = re.compile('[' + ''.join(dict_.keys()) + ']')
    num_string = patt.sub(lambda m: dict_[m.group(0)] + ' ', string)
    num = [int(i) for i in num_string.split()]
    return num

def read_record(file_, num_evo_entries):
    """ Read a Mathematica protein record from file and convert into dict. """
    
    dict_ = {}

    while True:
        next_line = file_.readline()
        for case in switch(next_line):
            if case('[ID]' + '\n'):
                id_ = file_.readline()[:-1]
                dict_.update({'id': id_})
            elif case('[PRIMARY]' + '\n'):
                primary = letter_to_num(file_.readline()[:-1], _aa_dict)
                dict_.update({'primary': primary})
            elif case('[EVOLUTIONARY]' + '\n'):
                pos = file_.tell()
                line = file_.readline()
                file_.seek(pos)
                evolutionary = []
                # Workaround for records with >20 rows in the PSSM
                while "[" not in line:
                    evolutionary.append([float(step) for step in file_.readline().split()])
                    pos = file_.tell()
                    line = file_.readline()
                    file_.seek(pos)
                dict_.update({'evolutionary': evolutionary})
            elif case('[SECONDARY]' + '\n'):
                secondary = letter_to_num(file_.readline()[:-1], _dssp_dict)
                dict_.update({'secondary': secondary})
            elif case('[TERTIARY]' + '\n'):
                tertiary = []
                for axis in range(NUM_DIMENSIONS): tertiary.append([float(coord) for coord in file_.readline().split()])
                dict_.update({'tertiary': tertiary})
            elif case('[MASK]' + '\n'):
                mask = letter_to_num(file_.readline()[:-1], _mask_dict)
                dict_.update({'mask': mask})
            elif case('[BFACTOR]' + '\n'):
                bfactor = [float(bf) for bf in file_.readline().split()]
                dict_.update({'bfactor': bfactor})
            elif case('\n'):
                return dict_
            elif case(''):
                return None

def dict_to_tfrecord(dict_):
    """ Convert protein dict into TFRecord. """

    id_ = _bytes_feature([dict_['id']])

    feature_lists_dict = {}
    feature_lists_dict.update({'primary': _feature_list([_int64_feature([aa]) for aa in dict_['primary']])})
    
    if dict_.has_key('evolutionary'): 
        feature_lists_dict.update({'evolutionary': _feature_list([_float_feature(list(step)) for step in zip(*dict_['evolutionary'])])})

    if dict_.has_key('secondary'):     
        feature_lists_dict.update({'secondary': _feature_list([_int64_feature([dssp]) for dssp in dict_['secondary']])})

    if dict_.has_key('tertiary'): 
        feature_lists_dict.update({'tertiary': _feature_list([_float_feature(list(coord)) for coord in zip(*dict_['tertiary'])])})
    
    if dict_.has_key('mask'):
        feature_lists_dict.update({'mask': _feature_list([_float_feature([step]) for step in dict_['mask']])})

    if dict_.has_key('bfactor'):
        feature_lists_dict.update({'bfactor': _feature_list([_float_feature([bf]) for bf in dict_['bfactor']])})


    record = _sequence_example(context=_features({'id': id_}), feature_lists=_feature_lists(feature_lists_dict))
    
    return record

# main. accepts one command-line argument: input file
# writes tf records of length 256 named 1...n
# activate tf-gpu environment
if __name__ == '__main__':
    input_path      = sys.argv[1] 
    output_path = ""
    num_evo_entries = 0

    input_file = open(input_path, 'r')

    protein_counter = 0
    output_counter = 1
    output_counter_string = str(output_counter)
    output_file = tf.io.TFRecordWriter(output_counter_string) 
    
    while True:
        dict_ = read_record(input_file, num_evo_entries)
        if dict_ is not None:
            tfrecord_serialized = dict_to_tfrecord(dict_).SerializeToString()
            print(protein_counter)
            protein_counter = protein_counter + 1
            output_file.write(tfrecord_serialized)
            if protein_counter % 256 == 0:
                output_counter = output_counter + 1
                output_counter_string = str(output_counter)
                output_file.close()
                output_file = tf.io.TFRecordWriter(output_counter_string)
        else:
            input_file.close()
            output_file.close()
            break
