# Beam Steering- Recursive Least Square

This is the [UMass Lowell SPACE HAUC](https://www.uml.edu/Research/LoCSST/Research/spacehauc/about.aspx) Beam Steering team source code.

#### Key Points:

* RLS is an optimization algorithm which recursively finds the coefficients to minimize the weighted cost function (essentially the finds coefficients that make the weights of the array more correct using previous data)

  _For example: x is the signal we receive, d is the desired signal, e is error, and w is the weighted coeff._
  
  ![image](https://user-images.githubusercontent.com/9031240/28492512-29762d08-6ed3-11e7-86c7-c677e4b1612c.png)

* RLS is a deterministic while the LMS is a stochastic function. Stochastic functions are randomly determined. Deterministic models involve no randomness in the ongoing development of the system, thus always produce the same output from a given starting condition. 
* RLS is more suited towards a fast changing mobile environment due to the slow speed of the LMS convergence. This problem is solved with replacing the gradient step size with a gain matrix at each iteration allowing it to update the weight estimates recursively as new data arrives. This makes the method more computationally efficient but more complex. 

##### The major pros of this method are that it is more computationally efficient, has a faster convergence speed, and is deterministic. It is suited for a fast, changing environment which is what we are dealing with. The cons are that it is more complex computationally and conceptually



###### This repository is filed under the GNU GPL v3 license in the spirit of promoting open source, and is free for anyone to use in their own projects.

[![GPL License](http://darrienglasser.com/gpl-v3-logo.jpg)](http://www.gnu.org/licenses/gpl-3.0.en.html)
