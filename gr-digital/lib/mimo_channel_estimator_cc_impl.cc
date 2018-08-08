/* -*- c++ -*- */
/* 
 * Copyright 2018 Free Software Foundation, Inc.
 * 
 * This file is part of GNU Radio
 * 
 * GNU Radio is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3, or (at your option)
 * any later version.
 * 
 * GNU Radio is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with GNU Radio; see the file COPYING.  If not, write to
 * the Free Software Foundation, Inc., 51 Franklin Street,
 * Boston, MA 02110-1301, USA.
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <gnuradio/io_signature.h>
#include "mimo_channel_estimator_cc_impl.h"
#include <pmt/pmt.h>
#define _USE_MATH_DEFINES
#include <math.h>

#if (EIGEN3_ENABLED)
#include <eigen3/Eigen/Dense>
#endif

#include <boost/format.hpp>
using namespace boost;

namespace gr {
  namespace digital {

    const pmt::pmt_t mimo_channel_estimator_cc_impl::d_key = pmt::string_to_symbol("pilot");

    mimo_channel_estimator_cc::sptr
    mimo_channel_estimator_cc::make(uint16_t M, uint16_t N, std::vector<std::vector<gr_complex> > training_sequence)
    {
      return gnuradio::get_initial_sptr
        (new mimo_channel_estimator_cc_impl(M, N, training_sequence));
    }

    /*
     * The private constructor
     */
    mimo_channel_estimator_cc_impl::mimo_channel_estimator_cc_impl(uint16_t M, uint16_t N, std::vector<std::vector<gr_complex> > training_sequence)
      : gr::block("mimo_channel_estimator_cc",
              gr::io_signature::make(N, N, sizeof(gr_complex)),
              gr::io_signature::make(N, N, sizeof(gr_complex))),
        d_M(M), d_N(N), d_training_sequence(training_sequence)
    {
      d_training_length = d_training_sequence[0].size();
      d_csi = std::vector<std::vector<gr_complex> >(N, std::vector<gr_complex> (M, 1.0));
      // Set the minimum size for the input buffer to the training length.
      set_min_noutput_items(d_training_length);
    }

    /*
     * Our virtual destructor.
     */
    mimo_channel_estimator_cc_impl::~mimo_channel_estimator_cc_impl()
    {
    }

    void
    mimo_channel_estimator_cc_impl::forecast (int noutput_items, gr_vector_int &ninput_items_required)
    {
      for (int i = 0; i < d_N; ++i) {
        ninput_items_required[i] = noutput_items;
      }
    }

    void
    mimo_channel_estimator_cc_impl::copy_symbols(gr_vector_const_void_star &input_items,
                                                 gr_vector_void_star &output_items,
                                                 uint32_t symbol_length,
                                                 uint32_t reading_offset,
                                                 uint32_t writing_offset) {
      for (int i = 0; i < d_N; ++i) {
        const gr_complex *in = (const gr_complex *) input_items[i];
        gr_complex *out = (gr_complex *) output_items[i];
        memcpy(&out[writing_offset], &in[reading_offset], symbol_length);
      }
    }

    void
    mimo_channel_estimator_cc_impl::estimate_channel(gr_vector_const_void_star &input_items,
                                                     uint32_t reading_offset) {
      switch (d_M) {
        case 1: {
          // Init CSI vector.
          std::vector<std::vector<gr_complex> > csi (d_N, std::vector<gr_complex> (1, 0.0));
          // Fill CSI vector with precalculated ML channel estimation.
          for (int n = 0; n < d_N; ++n) {
            const gr_complex* in = (const gr_complex *) input_items[n];
            for (int i = 0; i < d_training_length; ++i) {
              csi[n][0] += in[i];
            }
            // Normalize elements referring the training length.
            csi[n][0] *= 1./ d_training_length;
          }
          // Write local vector to class member.
          d_csi = csi;
          break;
        }
        case 2: {
          // Init csi vector.
          std::vector<std::vector<gr_complex> > csi (d_N, std::vector<gr_complex> (2, 0.0));
          // Fill CSI vector with precalculated ML channel estimation.
          for (int n = 0; n < d_N; ++n) {
            const gr_complex* in = (const gr_complex *) input_items[n];
            for (int i = 0; i < d_training_length; ++i) {
              csi[n][0] += in[i];
              csi[n][1] += in[i] * (gr_complex) std::polar(1.0, 2 * M_PI * i / d_training_length);
            }
            // Normalize elements referring the training length.
            csi[n][0] *= 1./ d_training_length;
            csi[n][1] *= 1./ d_training_length;
          }
          // Write local vector to class member.
          d_csi = csi;
          break;
        }
        default: {
#if (EIGEN3_ENABLED)
          // TODO calculate pseudo-inverse in constructor
          // Map training sequence and received sequence 2-dimensional std::vector to Eigen MatrixXcf.
          Eigen::MatrixXcf training_matrix(d_M, d_training_length);
          Eigen::MatrixXcf received_training_matrix(d_N, d_training_length);
          for (int i = 0; i < d_M; ++i) {
              training_matrix.row(i) = Eigen::VectorXcf::Map(&d_training_sequence[i][0], d_training_length);
          }
          for (int i = 0; i < d_N; ++i) {
              received_training_matrix.row(i) = Eigen::VectorXcf::Map((const gr_complex *) input_items[i], d_training_length);
          }
          // Calculate the Maximum Likelihood estimation for the channel matrix.
          Eigen::MatrixXcf csi = received_training_matrix * training_matrix.adjoint() * (training_matrix*training_matrix.adjoint()).inverse();

           // Map the CSI Eigen MatrixXcf to a 2-dim std::vector.
          for (int i = 0; i < d_N; ++i) {
            Eigen::VectorXcf::Map(&d_csi[i][0], d_M) = csi.row(i);
          }
#else
          throw std::runtime_error("Required library Eigen3 for MxN MIMO schemes with M,N>2 not installed.");
#endif
        }
      }
    }

    int
    mimo_channel_estimator_cc_impl::general_work (int noutput_items,
                       gr_vector_int &ninput_items,
                       gr_vector_const_void_star &input_items,
                       gr_vector_void_star &output_items)
    {
      // Collect all tags of the input buffer with key "csi" in the vector 'tags'.
      get_tags_in_window(tags, 0, 0, noutput_items, d_key);

      uint32_t nconsumed = 0; // Number of read items.
      uint32_t nwritten = 0; // Number of written items.
      uint32_t symbol_length; // Number of items in the current symbol.

      if(tags.size() == 0){ // Input buffer includes no tags at all.
        // Handle all samples in buffer as they belong to the current symbol.
        symbol_length = noutput_items;
        // Copy the data to the output buffer.
        copy_symbols(input_items, output_items, symbol_length, 0, 0);
        nconsumed = noutput_items;
        nwritten = noutput_items;
      } else { // Input buffer includes tags.
        if (tags[0].offset - nitems_read(0) > 0){
          /* There are items in the input buffer, before the first tag arrives,
+          * which belong to the previous symbol. */
          symbol_length = tags[0].offset - nitems_read(0);
          // Copy the data to the output buffer.
          copy_symbols(input_items, output_items, symbol_length, 0, 0);
          nconsumed += symbol_length;
          nwritten += symbol_length;
        }
        // Iterate over tags in buffer.
        for (unsigned int i = 0; i < tags.size(); ++i) {
          // Calculate the number of items before the next tag.
          if (i < tags.size() - 1) {
            // This is not the last tag in the buffer.
            symbol_length = tags[i + 1].offset - nitems_read(0) - nconsumed;
            // Check if the symbol length is smaller than the training length
            if (symbol_length < d_training_length){
              // We ignore this unfinished training sequence and dump it without any estimation.
              nconsumed += symbol_length;
              GR_LOG_INFO(d_logger, format("The symbol length (%d) is smaller than the training length (%d)."
                                           "No estimation possible. Dumping symbol.")
                                    %symbol_length %d_training_length);
              continue;
            }
          } else {
            // This is the last tag in the buffer.
            symbol_length = noutput_items - nconsumed;
            // Check if the training sequence of this last symbol is completely in this buffer.
            if (symbol_length < d_training_length){
              /* We consume the whole training sequence in the next buffer and don't consume it now
               * in order to do the estimation at once. */
              break;
            }
          }
          // Estimate MIMO channel and write CSI tag to stream.
          estimate_channel(input_items, nconsumed);

          // Consume training symbols. We don't copy them to the output buffer.
          nconsumed += d_training_length;
          symbol_length -= d_training_length;

          // Assign the CSI vector to a PMT vector.
          pmt::pmt_t csi_pmt = pmt::make_vector(d_N, pmt::make_c32vector(d_M, d_csi[0][0]));
          for (int n = 0; n < d_N; ++n){
            pmt::pmt_t csi_line_vector = pmt::make_c32vector(d_M, d_csi[n][0]);
            for (int m = 0; m < d_M; ++m) {
              pmt::c32vector_set(csi_line_vector, m, d_csi[n][m]);
            }
            pmt::vector_set(csi_pmt, n, csi_line_vector);
          }

          // Copy symbols to output and append stream tags with CSI to the beginning of this symbol.
          copy_symbols(input_items, output_items, symbol_length, nconsumed, nwritten);
          add_item_tag(0, nitems_written(0) + nwritten, pmt::string_to_symbol(std::string("csi")), csi_pmt);
          nconsumed += symbol_length;
          nwritten += symbol_length;
        }
      }
      // Tell runtime system how many input items we consumed on each input stream.
      consume_each (nconsumed);
      // Tell runtime system how many output items we produced.
      return nwritten;
    }

  } /* namespace digital */
} /* namespace gr */

