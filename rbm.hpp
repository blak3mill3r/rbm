/*
 * =====================================================================================
 *
 *       Filename:  rbm.hpp
 *
 *    Description   RBMPL main include file
 *
 *        Version:  1.0
 *        Created:  02/06/2009 03:00:47 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Blake Miller (bnm), blak3mill3r@gmail.com
 *        Company:  Stupendous Software
 *
 * =====================================================================================
 */

#include <iostream>
#include <math.h>

using std::cout;
using std::endl;

namespace rbm {

typedef float weight_t;
typedef unsigned int uint;

float sigmoid_immediate( float x )
{
  float a = pow(M_E, -x);         // M_E is e, the base of natural logarithms
  return ( 1 / ( 1 + a ) );       // the function is f(x) = 1 / (1 + e^-x)
}

float sigmoid( float x )
{ return sigmoid_immediate( x ); }

template <class ActivationType>
class Neuron
{
  public:
    typedef ActivationType activation_type;
    ActivationType    activation;
    Neuron()
    { activation = 0; }

    inline void compute_activation( weight_t a )
    {
      // stochastic binary neuron activation from sigmoid function
    }
};

template <>
Neuron< bool >::Neuron()
{ activation = false; }

template <>
Neuron< float >::Neuron()
{ activation = 0.1; }

template < class NeuronType >
class Neurons    // a blob of neurons or perceptrons
{
  public:
    NeuronType        *neurons;

    uint              num_neurons;

    typedef Neurons< NeuronType > self;
    typedef NeuronType neuron_type;

    self * connectable_neurons()
    { return this; }
    
    uint num_connectable_neurons()
    { return num_neurons; }

    Neurons( uint n ) : num_neurons(n)
    {
      cout << "Neurons<class NeuronType>::Neurons("<< n << ")" << endl;
      neurons = new NeuronType[ num_neurons ];
    }

    virtual ~Neurons()
    { delete [] neurons; }

};

// handles 1 layer of RBMness
template <
    class Lower    // type of lower-level Neurons or RBM
    , class Upper  // type of upper-level Neurons
>
class RBM       // Restricted Boltzman Majigger
{
  public:
    //typedef RBM< Lower, Upper > self;
    typedef typename Upper::neuron_type neuron_type;

    Lower    *lower;
    Upper    *upper;

    uint              num_lower_neurons, num_upper_neurons;

    Upper * connectable_neurons()
    { return upper; }

    uint num_connectable_neurons()
    { return num_upper_neurons; }

    RBM( Lower * lower_, uint num_upper_neurons_ )
    {
      cout << "RBM()" << endl;
      // allocate memory for weights
      // allocate an upper layer
      lower = lower_;
      num_lower_neurons = lower->num_connectable_neurons();
      cout << "weights =new weight_t[" << (num_lower_neurons * num_upper_neurons_ ) << "];" << endl;
      weights = new weight_t[ num_lower_neurons * num_upper_neurons_ ];
      num_upper_neurons = num_upper_neurons_;
      upper = new Upper( num_upper_neurons );

      // FIXME optionally allocate lower layer or lower RBM?
    }
    
    virtual ~RBM()
    {
      cout << "~RBM()" << endl;
      // deallocate weights & upper layer
      delete [] weights;
      delete upper;
    }

    /*void train( const what* data, int n_frames, const what* labels = NULL )
    {}

    void perceive( const what* data )
    {}

    void fantasize( const what* labels = NULL )
    {}
*/


  private:
    weight_t *weights;

    // for one neuron in our upper level
    // activate according to its inputs from the lower level
    inline void upward_activation_step( uint upper_neuron_id )
    {
      static weight_t sum = 0;
      uint lower_neuron_id = 0, weight_id = upper_neuron_id * num_lower_neurons;

      while(lower_neuron_id < num_lower_neurons)
        sum += lower->neurons[lower_neuron_id++] * weights[weight_id++];
      
      upper->neurons[ upper_neuron_id ].activation = sigmoid( sum );
    }

    // for one neuron in our lower level
    // activate according to its outputs from the higher level
    inline void downward_activation_step( uint lower_neuron_id )
    {
      static weight_t sum = 0;
      uint upper_neuron_id = 0, weight_id = lower_neuron_id;

      while(upper_neuron_id < num_upper_neurons)
      {
        sum += upper->neurons[upper_neuron_id++] * weights[weight_id];
        weight_id += num_lower_neurons;
      }

      lower->connectable_neurons()->neurons[ lower_neuron_id ].activation = sigmoid( sum );
    }

};

} // end namespace rbm
