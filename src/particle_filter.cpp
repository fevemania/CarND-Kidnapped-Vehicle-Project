/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  // TODO: Set the number of particles. Initialize all particles to first position (based on estimates of
  //   x, y, theta and their uncertainties from GPS) and all weights to 1.
  // Add random Gaussian noise to each particle.
  // NOTE: Consult particle_filter.h for more information about this method (and others in this file).
  num_particles = 50;

  normal_distribution<double> noise_x(x, std[0]);
  normal_distribution<double> noise_y(y, std[1]);
  normal_distribution<double> noise_theta(theta, std[2]);

  for (int i = 0; i < num_particles; ++i) {
    Particle tmp;
    tmp.x = noise_x(dre);
    tmp.y = noise_y(dre);
    tmp.theta = noise_theta(dre);
    weights.push_back(1.0);

    particles.push_back(std::move(tmp));
  }

  is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
  // TODO: Add measurements to each particle and add random Gaussian noise.
  // NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
  //  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
  //  http://www.cplusplus.com/reference/random/default_random_engine/

  normal_distribution<double> noise_x(0.0, std_pos[0]);
  normal_distribution<double> noise_y(0.0, std_pos[1]);
  normal_distribution<double> noise_theta(0.0, std_pos[2]);

  for (auto& particle : particles) {
    double pred_theta = particle.theta + yaw_rate * delta_t;
    particle.x += (velocity / yaw_rate) * ( sin(pred_theta) - sin(particle.theta) ) + noise_x(dre);
    particle.y += (velocity / yaw_rate) * ( cos(particle.theta) - cos(pred_theta) ) + noise_y(dre);
    particle.theta = pred_theta + noise_theta(dre);
  }
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
  // TODO: Find the predicted measurement that is closest to each observed measurement and assign the
  //   observed measurement to this particular landmark.
  // NOTE: this method will NOT be called by the grading code. But you will probably find it useful to
  //   implement this method and use it as a helper during the updateWeights phase.


}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
                                   const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
  // TODO: Update the weights of each particle using a mult-variate Gaussian distribution with respect to observations.
  //   You can read more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
  // NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
  //   according to the MAP'S coordinate system. You will need to transform between the two systems.
  //   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
  //   The following is a good resource for the theory:
  //   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
  //   and the following is a good resource for the actual equation to implement (look at equation
  //   3.33
  //   http://planning.cs.uiuc.edu/node99.html

  // for each particle

  for (int i = 0; i < particles.size(); ++i) {

    // initial weight
    double weight = 1.0;

    for (int j = 0; j < observations.size(); ++j) {
      // Step1. translate car sensor landmark obs from VEHICLE's coordinate system to MAP's coordinate system
      LandmarkObs ob_in_map;
      ob_in_map.x = particles[i].x +
              (cos(particles[i].theta) * observations[j].x - sin(particles[i].theta) * observations[j].y);
      ob_in_map.y = particles[i].y +
              (sin(particles[i].theta) * observations[j].x + cos(particles[i].theta) * observations[j].y);

      // Step2. associate with the closest map landmark
      double min_distance = numeric_limits<double>::max();
      for (const auto& mark: map_landmarks.landmark_list) {
        double dist_tmp = dist(ob_in_map.x, ob_in_map.y, mark.x_f, mark.y_f);
        if (dist_tmp < sensor_range && dist_tmp < min_distance) {
          min_distance = dist_tmp;
          //  ob_in_map.id = mark.id_i;
        }
      }

      // Step3. calculate particle's weight
      double var_x = pow(std_landmark[0], 2);
      double var_y = pow(std_landmark[1], 2);
      double covar_xy = std_landmark[0] * std_landmark[1];

      double exponent = (min_distance * min_distance / 2 * var_x * var_y);

      // final weight is the multiplication of all calculated measurement probabilities from observations.
      weight *= (exp(-exponent) / 2*M_PI * covar_xy);
    }
    particles[i].weight = weights[i] = weight;
  }


}

void ParticleFilter::resample() {
  // TODO: Resample particles with replacement with probability proportional to their weight.
  // NOTE: You may find std::discrete_distribution helpful here.
  //   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
  std::vector<Particle> newParticles;
  std::discrete_distribution<int> sample_according_to_weight(weights.begin(), weights.end()); // return index to particles

  for (const auto& _ : particles)
    newParticles.push_back(particles[sample_according_to_weight(dre)]);

  particles = std::move(newParticles);
}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations,
                                         const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
  // particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best)
{
  vector<int> v = best.associations;
  stringstream ss;
  copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}
string ParticleFilter::getSenseX(Particle best)
{
  vector<double> v = best.sense_x;
  stringstream ss;
  copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}
string ParticleFilter::getSenseY(Particle best)
{
  vector<double> v = best.sense_y;
  stringstream ss;
  copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}
