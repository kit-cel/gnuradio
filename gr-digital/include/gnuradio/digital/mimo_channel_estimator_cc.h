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


#ifndef INCLUDED_DIGITAL_MIMO_CHANNEL_ESTIMATOR_CC_H
#define INCLUDED_DIGITAL_MIMO_CHANNEL_ESTIMATOR_CC_H

#include <gnuradio/digital/api.h>
#include <gnuradio/block.h>

namespace gr {
  namespace digital {

    /*! \brief Estimates a MIMO channel matrix.
     *
     * The block estimates the channel matrix of a MIMO scheme with the help of
     * training sequences. The training sequence for each transmitting antenna
     * equals one subvector/row of the 2-dimensional vector training_sequence.
     * The training sequence should be appended as a pilot to the data at the transmitter.
     * The beginning of the pilot must be tagged at the receiver with the key 'pilot'.
     * This tag can be set from a preceding sync block, for example.
     *
     * This block has N inports and N outports. It estimates the channel matrix which each incoming
     * 'pilot' tag, dumps the pilot symbol of each stream and passes the rest of the data through without
     * changing anything. A tag is set at the beginning of each symbol (= a data sequence between to training sequences)
     * with the key 'csi' that contains the estimated MxN channel matrix.
     *
     * For 1xN and 2xN MIMO schemes, the equalizer is calculated internally. For MxN schemes with M > 2,
     * the C++ template library Eigen is used for the linear algebra operations
     * and is a requirement in this case.
     * The training_length is not bounded above, but the minimum length is M
     * (to avoid an underdetermined equation system).
     */
    class DIGITAL_API mimo_channel_estimator_cc : virtual public gr::block
    {
     public:
      typedef boost::shared_ptr<mimo_channel_estimator_cc> sptr;

      /*!
       * \brief Return a shared_ptr to a new instance of digital::mimo_channel_estimator_cc.
       *
       * To avoid accidental use of raw pointers, digital::mimo_channel_estimator_cc's
       * constructor is in a private implementation
       * class. digital::mimo_channel_estimator_cc::make is the public interface for
       * creating new instances.
       */
      static sptr make(uint16_t M, uint16_t N, std::vector<std::vector<gr_complex> > training_sequence);
    };

  } // namespace digital
} // namespace gr

#endif /* INCLUDED_DIGITAL_MIMO_CHANNEL_ESTIMATOR_CC_H */

