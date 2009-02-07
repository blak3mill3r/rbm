/*
 * =====================================================================================
 *
 *       Filename:  rbmpl_test.cpp
 *
 *    Description:  simple test for RBMPL
 *
 *        Version:  0.01
 *        Created:  02/06/2009 02:58:13 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Blake Miller (bnm), blak3mill3r@gmail.com
 *        Company:  Stupendous Software
 *
 * =====================================================================================
 */

#include <jack/jack.h>
#include <iostream>
#include "rbmpl.hpp"

using namespace rbmpl;

using std::cout;
using std::endl;

jack_port_t *input_port;
jack_client_t *client;
jack_default_audio_sample_t *samples_of_a = NULL;

  typedef Neuron< jack_default_audio_sample_t > audio_neuron_type;
  typedef Neurons< audio_neuron_type > audio_neurons_type;

  typedef Neuron< bool > neuron_type;
  typedef Neurons< neuron_type > neurons_type;

  typedef RBM< audio_neurons_type, neurons_type > first_level_rbm_type;
  typedef RBM< first_level_rbm_type,  neurons_type > second_level_rbm_type;

int main(int argc, char ** argv) {
  audio_neurons_type audio_perceptrons( 512 );
  cout  << "made audio_perceptrons at " << &audio_perceptrons << endl;
  cout << "audio_perceptrons.num_connectable_neurons() = " << audio_perceptrons.num_connectable_neurons() << endl;
  first_level_rbm_type sub_net( &audio_perceptrons, 512 );
  cout  << "made sub_net at " << &sub_net << endl;
  cout << "sub_net.num_connectable_neurons() = " << sub_net.num_connectable_neurons() << endl;
  second_level_rbm_type net( &sub_net, 512 );
  cout  << "made net at " << &net << endl;
  cout << "net.num_connectable_neurons() = " << net.num_connectable_neurons() << endl;
}

