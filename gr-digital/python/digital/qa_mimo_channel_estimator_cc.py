#!/usr/bin/env python
# -*- coding: utf-8 -*-
# 
# Copyright 2018 Free Software Foundation, Inc.
# 
# This file is part of GNU Radio
# 
# GNU Radio is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3, or (at your option)
# any later version.
# 
# GNU Radio is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with GNU Radio; see the file COPYING.  If not, write to
# the Free Software Foundation, Inc., 51 Franklin Street,
# Boston, MA 02110-1301, USA.
# 

from gnuradio import gr, gr_unittest
from gnuradio import blocks
import digital_swig as digital
import numpy as np
import pmt

class qa_mimo_channel_estimator_cc (gr_unittest.TestCase):

    def setUp (self):
        self.tb = gr.top_block ()

    def tearDown (self):
        self.tb = None

    def produce_training_seq(self, M, length):
        seq = np.empty([M, length], dtype=complex)
        for (m, i), value in np.ndenumerate(seq):
            seq[m][i] = np.exp(-1j*2*np.pi*m*i/length)
        return seq

    def run_flowgraph(self, M, N, training_sequence, training_length, channel_matrix):
        # Produce tags.
        tags = [gr.tag_utils.python_to_tag((0,
                                            pmt.string_to_symbol("pilot"),
                                            pmt.from_long(0),
                                            pmt.from_long(0)))]
        # Build GNU Radio blocks.
        data = np.random.randn(training_length) + 1j * np.random.randn(training_length)
        src = [blocks.vector_source_c(data=np.append(training_sequence[0], data),
                                      repeat=False,
                                      tags=tags)]
        sink = [blocks.vector_sink_c_make()]
        for m in range(1, M):
            src.append(blocks.vector_source_c(np.append(training_sequence[m], data)))
        for n in range(1, N):
            sink.append(blocks.null_sink_make(gr.sizeof_gr_complex))
        # Use a matrix multiplier as a static channel model.
        channel = blocks.multiply_matrix_cc_make(channel_matrix)
        estimator = digital.mimo_channel_estimator_cc(M, N, training_sequence)
        # Use a head block to terminate the flowgraph.
        head = blocks.head_make(gr.sizeof_gr_complex, training_length)
        # Connect everything.
        for m in range(0, M):
            self.tb.connect(src[m], (channel, m))
        self.tb.connect((channel, 0), (estimator, 0), head, sink[0])
        for n in range(1, N):
            self.tb.connect((channel, n), (estimator, n), sink[n])
        # Run flowgraph.
        self.tb.run()
        # Read channel estimation and write from pmt to numpy array.
        csi = np.empty(shape=[N, M], dtype=complex)
        for n in range(0, N):
            for m in range(0, M):
                csi[n][m] = pmt.c32vector_ref(pmt.vector_ref(sink[0].tags()[0].value, n), m)
        return csi


    '''
    2 tests validating the correct estimation of the channel coefficients with a random channel
    and 1xN MIMO scheme. '''
    def test_1xN_scheme (self):
        # Test parameters.
        M = 1
        repetitions = 2
        for i in range(0, repetitions):
            N = np.random.randint(1, 65)
            training_length = np.random.randint(M, 65)
            # Produce training sequence and random channel coefficients.
            training_sequence = self.produce_training_seq(M, training_length)
            channel_matrix = (np.random.randn(N, M) + 1j * np.random.randn(N, M))

            estimation = self.run_flowgraph(M, N, training_sequence, training_length, channel_matrix)

            # Check if the expected result equals the actual result.
            for n in range(0, N):
                self.assertComplexTuplesAlmostEqual(channel_matrix[n], estimation[n], 2)

    '''
    2 tests validating the correct estimation of the channel coefficients with a random channel
    and 2xN MIMO scheme. '''
    def test_2xN_scheme (self):
        # Test parameters.
        M = 2
        repetitions = 2
        for i in range(0, repetitions):
            N = np.random.randint(1, 65)
            training_length = np.random.randint(M, 65)
            # Produce training sequence and random channel coefficients.
            training_sequence = self.produce_training_seq(M, training_length)
            channel_matrix = (np.random.randn(N, M) + 1j * np.random.randn(N, M))

            estimation = self.run_flowgraph(M, N, training_sequence, training_length, channel_matrix)

            # Check if the expected result equals the actual result.
            for n in range(0, N):
                self.assertComplexTuplesAlmostEqual(channel_matrix[n], estimation[n], 2)

    '''
    2 tests validating the correct estimation of the channel coefficients with a random channel
    and MxN MIMO scheme (M>2). '''
    def test_MxN_scheme (self):
        # Test parameters.
        M = np.random.randint(3, 65)
        repetitions = 2
        for i in range(0, repetitions):
            N = np.random.randint(1, 65)
            training_length = np.random.randint(M, 65)
            # Produce training sequence and random channel coefficients.
            training_sequence = self.produce_training_seq(M, training_length)
            channel_matrix = (np.random.randn(N, M) + 1j * np.random.randn(N, M))

            estimation = self.run_flowgraph(M, N, training_sequence, training_length, channel_matrix)

            # Check if the expected result equals the actual result.
            for n in range(0, N):
                self.assertComplexTuplesAlmostEqual(channel_matrix[n], estimation[n], 2)

if __name__ == '__main__':
    gr_unittest.run(qa_mimo_channel_estimator_cc, "qa_mimo_channel_estimator_cc.xml")
