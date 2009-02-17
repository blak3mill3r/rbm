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
#include <iomanip>
#include <math.h>

using std::cout;
using std::cerr;
using std::endl;
using std::setw;
using std::left;
using std::right;
using std::showpos;
using std::noshowpos;

namespace rbm {

typedef float weight_t;
typedef unsigned int uint;

float sigmoid_immediate( float x )
{
  float temperature = 1.0;
  float a = pow(M_E, (-x * temperature));         // M_E is e, the base of natural logarithms
  return ( 1 / ( 1 + a ) );       // the function is f(x) = 1 / (1 + e^-x)
}

float sigmoid( float x )          // FIXME optimize
{ return sigmoid_immediate( x ); }

template <class ActivationType>
class Neuron
{
  public:
    typedef ActivationType activation_type;
    weight_t          bias;

    Neuron()
      { m_activation = 0; }

    inline weight_t activation()
      { return m_activation; }

    inline void set_activation( activation_type a )
      { m_activation = a; }

  private:
    ActivationType    m_activation;
};

template <>
Neuron< float >::Neuron()
{ m_activation = 0.1; }

template <class ActivationType>
class StochasticNeuron              // fires or doesn't fire, like a real neuron
{
  public:
    typedef ActivationType activation_type;
    weight_t          bias;

    StochasticNeuron()
      { m_activation = 0; }

    inline weight_t activation()
      {
        //decide!
        float r = ((float)rand())/RAND_MAX;
        //cout << m_activation << '*' << r << endl;
        if ( r > ((m_activation + 1.0)*0.5) ) 
          return 1.0;
        else
          return 0.0;
      }

    inline void set_activation( activation_type a )
      { m_activation = a; }

  private:
    ActivationType    m_activation;
};


// primary template for class Neurons
// description: a blob of neurons or perceptrons
//   passing something to the second argument of the constructor will cause it to be a "perceptron" blob
//   constructed with one argument it's just a blob of neurons of a given NeuronType, the resources for neuron values are managed by this class
template < class NeuronType >
class Neurons
{
  public:
    NeuronType        *neurons;

    uint              num_neurons;

    typedef Neurons< NeuronType > self;
    typedef NeuronType neuron_type;
    typedef typename NeuronType::activation_type activation_type;

///FIXME are these next two function necessary??
    self * connectable_neurons()
    { return this; }
    
    uint num_connectable_neurons()
    { return num_neurons; }

    void set_activation( activation_type * v )
    {
      for(uint i = 0; i < num_neurons; ++i)
        neurons[i].set_activation(v[i]);
    }

    Neurons( uint n ) : num_neurons(n)
    {
      //cout << "Neurons<class NeuronType>::Neurons("<< n << ")" << endl;
      neurons = new NeuronType[ num_neurons ];
    }

    virtual ~Neurons()
    { 
      delete [] neurons; 
    }

    void debugify()
    {
      cout << "Neurons::debugify()" << endl; 
      for(uint i = 0; i < num_neurons; ++i)
        cout << i << ':' << neurons[i].activation()<< ' '<< endl;
      cout << endl;
    }

};

// primary template for class RBM
// description: handles 1 layer of RBMness
template <
    class Lower    // type of lower-level Neurons
    , class Upper  // type of upper-level Neurons
>
class RBM       // Restricted Boltzman Majigger
{
  public:
    //typedef RBM< Lower, Upper > self;
    typedef typename Upper::neuron_type neuron_type;

    Lower    *lower;
    Upper    *upper;

    uint              num_lower_neurons, num_upper_neurons, num_connections, num_weights;

    Upper * connectable_neurons()
    { return upper; }

    uint num_connectable_neurons()
    { return num_upper_neurons; }

    RBM( Lower * lower_, uint num_upper_neurons_, weight_t learning_rate_ )
    {
      // initialize learning rate
      learning_rate = learning_rate_;
      // allocate memory for weights
      lower = lower_;
      num_lower_neurons = lower->num_connectable_neurons();
      num_upper_neurons = num_upper_neurons_;
      num_weights = num_lower_neurons * num_upper_neurons_ ;
      weights = new weight_t[ num_weights ];
      // allocate an upper layer
      upper = new Upper( num_upper_neurons );
    }
    
    virtual ~RBM()
    {
      delete [] weights;
      delete upper;
    }

    void randomize_weights()
    {
      srand(23897592);
      for(uint wi = 0; wi < num_weights; ++wi)
      {
        weights[wi] = (((float)rand())/RAND_MAX);
      }
    }

    void perceptual_learning(uint iterations)
    {
      for(uint i = 0; i < iterations; ++i)
        perceptual_learning_step();
    }

    // compare real data with reconstruction
    // this function assumes that the lower neurons have been set to real data
    void debugify()
    {
      // print real data
      cout << "the data:(skip)" << endl;
      //lower->debugify();
      // interpret it
      cout << "upact:" << endl;
      upward_activation();
      upper->debugify();
      // fantasize about it
      cout << "downact:" << endl;
      downward_activation();
      // print fantasy
      cout << "model's fantasy data:" << endl;
      lower->debugify();
    }

  private:
    weight_t *weights;
    weight_t learning_rate;

    // this function assumes that the lower neurons have been set to real data
    // it adjusts the weights by sampling from the net's energy distribution
    // in 2 states: after "interpreting" the data (upward activation step)
    // and after interpreting it's (fantasy based on it's interpretation of the data)
    // (up-down-up activation)
    void perceptual_learning_step()
    {
      upward_activation();

      // positive weight learning step:
      uint weight_id = 0;
      for( uint ui = 0; ui < num_upper_neurons; ++ui )
      {
        for( uint li = 0; li < num_lower_neurons; ++li )
        {
          // increment the weight by how active both connected neurons are
          weights[weight_id++] += (lower->neurons[li].activation() *
                                   upper->neurons[ui].activation() * 
                                   learning_rate);
        }
      }

      // positive bias learning step:
      /*for( uint li = 0; li < num_lower_neurons; ++li )
      {
        // increment the bias for the lower neuron
        lower->neurons[li].bias += lower->neurons[li].activation * learning_rate;
      }*/

      downward_activation();
      upward_activation();

      // negative weight learning step:
      weight_id = 0;
      for( uint ui = 0; ui < num_upper_neurons; ++ui )
      {
        for( uint li = 0; li < num_lower_neurons; ++li )
        {
          // decrement the weight by how active both connected neurons are
          weights[weight_id++] -= (lower->neurons[li].activation() *
                                   upper->neurons[ui].activation() * 
                                   learning_rate);
        }
      }

      // negative bias learning step:
      /*for( uint li = 0; li < num_lower_neurons; ++li )
      {
        // increment the bias for the lower neuron
        lower->neurons[li].bias -= lower->neurons[li].activation() * learning_rate;
      }*/

    } // end function alternating_gibbs_sample_step

    // for each neuron in our upper level
    // activate according to its inputs from the lower level
    void upward_activation()
    {
      for( uint upper_neuron_id = 0;
           upper_neuron_id < num_upper_neurons;
           ++upper_neuron_id )
        {
          weight_t sum = 0;
          uint lower_neuron_id = 0, weight_id = upper_neuron_id * num_lower_neurons;

          while(lower_neuron_id < num_lower_neurons)
            sum += lower->neurons[lower_neuron_id++].activation() * weights[weight_id++];
          
          upper->neurons[ upper_neuron_id ].set_activation( sigmoid(sum) );
        }
    }

    // for each neuron in our lower level
    // activate according to its inputs from the upper level
    void downward_activation()
    {
      //cout << "downward_activation()" << endl;
      for( uint lower_neuron_id = 0;
           lower_neuron_id < num_lower_neurons;
           ++lower_neuron_id )                                                                        // for each neuron in the lower set
        {
          //cout << "lower neuron[" << lower_neuron_id << ']';
          weight_t sum = 0.0f;                                                                    // sum the weights*activations of connected neurons in the upper set
          uint upper_neuron_id = 0, weight_id = lower_neuron_id;

          while(upper_neuron_id < num_upper_neurons)
          {
            //cout << "u[" << upper_neuron_id << ']';
            sum += upper->neurons[upper_neuron_id++].activation() * weights[weight_id];
            weight_id += num_lower_neurons;
            //cout << '=' << sum << '|' ;
          }

          lower->neurons[ lower_neuron_id ].set_activation( sigmoid(sum) );
        }

    }

};

} // end namespace rbm
