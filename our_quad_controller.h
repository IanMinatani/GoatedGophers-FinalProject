#include <cmath>

namespace aa448 {
	class quad_controller {
		private:
			float ex_num_array_1_[5][3]; // example data member that is a 5x3 array of floats (you can delete this line in your implementation).
			float ex_den_array_1_[5][3]; // example data member that is a 5x3 array of floats (you can delete this line in your implementation).
			// you should declare more internal data members here...
			// ...
			// ...
			// ...
			
		public:
			quad_controller() {
				// Constructor function that is called when quad_controller object
				// is instantiated.
				initialize(); // you should not need to change this.
			}
			~quad_controller() {
				// No-op. You should not need to change this.
			}
		
		public:
			void initialize() {
				// Initialize your data members.

				// Example (delete this in your implementation, and replace it with your own).
				memset(ex_num_array_1_,0,sizeof(float)*5*3); // memset sets the data pointed to by ex_num_array_1_,
				                                             // which has sizeof(float)*5*3 bytes of data, to zero. Without doing
				                                             // this, you are not guaranteed that ex_num_array_1_ will
				                                             // have zero values.

				// Example - an alternative to memset (delete this in your implementation).
				for (int i=0; i<5; ++i) {
					for (int j=0; j<3; ++j) {
						ex_den_array_1_[i][j] = 0;
					}
				}

				// You should initialize all of your data members here...
			}
			void step_wrench(const float r_f_z, const float r_torques_xyz[3], float pwms[4]) {
				// Interface:
				// 1. r_f_z         - 1x1 Input  - Reference z-force [N].
				// 2. r_torques_xyz - 3x1 Input  - Reference xyz-torque vector [N-m].
				// 3. pwms          - 4x1 Output - PWM commands [us].

				// This function runs the lowest-level part of your controller, and should contain
				// only the Control Allocator and the Signal Converter. It takes in the commanded
				// z-force and xyz-moments and converts them into PWM signals.
				//
				// This function is implemented for you, and you should not have to
				// modify it. However, you do need to implement the control_allocator()
				// and signal_converter() functions included further down.
				// 
				// This function has the interface indicated by the red dots in layout figure
				// in the project document.

				// Run control allocator.
				const float w_req[4] = {r_f_z,r_torques_xyz[0],r_torques_xyz[1],r_torques_xyz[2]};
				float f_mot_cmd[4] = {0,};
				control_allocator(w_req,f_mot_cmd);

				// Run signal converter.
				signal_converter(f_mot_cmd,pwms);
			}
			void step_rates(const float r_f_z, const float r_omega_xyz[3], const float omega_xyz_est[3], float pwms[4]) {
				// Interface:
				// 1. r_f_z         - 1x1 Input  - Reference z-force [N].
				// 2. r_omega_xyz   - 3x1 Input  - Reference angular velocity vector [rad/s].
				// 3. omega_xyz_est - 3x1 Input  - Measured angular velocity vector [rad/s].
				// 4. pwms          - 4x1 Output - PWM commands [us].
				//
				// This function runs your step_wrench() controller with an angular rate controller
				// wrapped around it. In this function, you should implement the C_omega_x/y/z,
				// F_omega_x/y/z, and P_omega_x/y/z blocks.
				//
				// You should leverage the flight controller you developed for Lab 5.
				// I would recommend starting simple, e.g., with a PID design. But
				// you can write your code in such a way that you can easily
				// move to a higher-order contorller. Reference the pid.h file in the updated Lab 6 code
				// for a function that converts kp, ki, kd, derivative LPF cutoff frequency,
				// and sampling time into digifilt numerator/denominator format.
				// 
				// The Controller Allocator and Signal Converter
				// blocks will be implemented implicitly by the call to the step_wrench()
				// function.
				//
				// This function has the interface indicated by the blue dots in layout figure
				// in the project document.

				// THIS IS A PLACEHOLDER TO PREVENT THE COMPILER FROM COMPLAINING.
				// DELETE THIS CODE IN YOUR IMPELEMENTATION AND REPLACE IT WITH
				// YOURS.
				(void)r_omega_xyz;
				(void)omega_xyz_est;
				const float r_torques_xyz[3] = {0,};

				// Call step_wrench().
				step_wrench(r_f_z,r_torques_xyz,pwms);
			}
			void step_attitude(const float r_f_z, const float q_i2r[4], const float omega_xyz_est[3], const float q_i2b_est[4], float pwms[4]) {
				// Interface:
				// 1. r_f_z         - 1x1 Input  - Reference z-force [N].
				// 2. q_i2r         - 4x1 Input  - Unit quaternion (wxyz) specifying the inertial-to-reference attitude.
				// 3. omega_xyz_est - 3x1 Input  - Measured angular velocity vector [rad/s].
				// 4. q_i2b_est     - 4x1 Input  - Unit quaternion (wxyz) specifying the inertial-to-body measured attitude.
				// 5. pwms          - 4x1 Output - PWM commands [us].
				//
				// This function runs your step_rates() controller with an attitude controller
				// wrapped around it. In this function, you should implement a proportional quaternion
				// controller that takes q_i2r and q_i2b_est and computes the reference angular velocity
				// r_omega_xyz.
				//
				// The angular-rate controller, Controller Allocator, and Signal Converter
				// blocks will be implemented implicitly by the call to the step_rates()
				// function.
				//
				// This function has the interface indicated by the green dots in layout figure
				// in the project document.

				// THIS IS A PLACEHOLDER TO PREVENT THE COMPILER FROM COMPLAINING.
				// DELETE THIS CODE IN YOUR IMPELEMENTATION AND REPLACE IT WITH
				// YOURS.
				(void)r_f_z;
				(void)q_i2r;
				(void)q_i2b_est;
				const float r_omega_xyz[3] = {0,};

				// Call step_rates().
				step_rates(r_f_z,r_omega_xyz,omega_xyz_est,pwms);
			}
			void step_4dof(const float r_a_z, const float q_i2r[4], const float omega_xyz_est[3], const float q_i2b_est[4], const float a_z_est, float pwms[4]) {
				// Interface:
				// 1. r_a_z         - 1x1 Input  - Reference z-acceleration [m/s^2].
				// 2. q_i2r         - 4x1 Input  - Unit quaternion (wxyz) specifying the inertial-to-reference attitude.
				// 3. omega_xyz_est - 3x1 Input  - Measured angular velocity vector [rad/s].
				// 4. q_i2b_est     - 4x1 Input  - Unit quaternion (wxyz) specifying the inertial-to-body measured attitude.
				// 5. a_z_est       - 1x1 Input  - Measured z-acceleration [m/s^2].
				// 6. pwms          - 4x1 Output - PWM commands [us].
				//
				// This function runs your step_attitude() controller with a z-acceleration controller
				// wrapped around it. In this function, you should implement a z-acceleration
				// feedback controller consisting of C_a_z, F_a_z, and P_a_z using code similar
				// to the one you employed in step_rates().
				//
				// The attitude controller, angular-rate controller, Controller Allocator,
				// and Signal Converter blocks will be implemented implicitly by the call
				// to the step_attitude() function.
				//
				// This function has the interface indicated by the orange dots in layout figure
				// in the project document.

				// THIS IS A PLACEHOLDER TO PREVENT THE COMPILER FROM COMPLAINING.
				// DELETE THIS CODE IN YOUR IMPELEMENTATION AND REPLACE IT WITH
				// YOURS.
				(void)r_a_z;
				(void)a_z_est;
				const float r_f_z = 0;

				// Call step_attitude().
				step_attitude(r_f_z,q_i2r,omega_xyz_est,q_i2b_est,pwms);
			}

		private:
			void control_allocator(const float w_req[4], float f_mot_cmd[4]) {
				// Interface:
				// 1. w_req     - Input  - 4x1 vector of wrench request: f_z [N], m_x [N-m], m_y [N-m], m_z [N-m].
				// 2. f_mot_cmd - Output - 4x1 vector of commanded motor forces: f_mot_cmd_1 [N], f_mot_cmd_2 [N], f_mot_cmd_3 [N], f_mot_cmd_4 [N].
				//
				// This function implements the control allocator. At a minimum, it
				// should invert the relationship w_req = A * f_mot_cmd, where A
				// is the 4x4 matrix that converts motor forces to total force and
				// xyz moments. That is, at a minimum, this function should perform
				// the operation f_mot_cmd = inv(A) * w_req. Note however that since
				// A is constant and known in advance, inv(A) can be computed offline.
				//
				// If you so desire, you can try to implement a more sophisticated
				// control allocator that prioritizes m_x and m_y as the top priority,
				// f_z as the second priority, and m_z as the lowest priority.

				float matrixHinverse[4][4] = { {-0.25, -2.9412, 2.9412, 16.6667},
				                               {-0.25, -2.9412, -2.9412, -16.6667}, 
    				               	           {-0.25, 2.9412, -2.9412, 16.6667},
											   {-0.25, 2.9412, 2.9412, -16.6667}};
				
				for (int i = 0; i < 4; i++) {
					f_mot_cmd[i] = 0;
				        for (int j = 0; j < 4; j++) {
				            f_mot_cmd[i] += matrixHinverse[i][j] * w_req[j];
				        }
				    }
				
			}
			void signal_converter(const float f_mot_cmd[4], float pwms[4]) {
				// Interface:
		                // 1. f_mot_cmd - Input  - 4x1 vector of commanded motor forces: f_mot_cmd_1 [N], f_mot_cmd_2 [N], f_mot_cmd_3 [N], f_mot_cmd_4 [N].
		                // 2. pwms      - Output - 4x1 vector of PWM signals: pwm_1 [us], pwm_2 [us], pwm_3 [us], pwm_4 [us]
		                float f[] = {f_mot_cmd[0],f_mot_cmd[1],f_mot_cmd[2],f_mot_cmd[3]};
		                float a = 1.5685e-6;
		                float b = 0.001126;
		                float c = 0.09276;
						const float f_min = 0.0381;
						const float f_max = 1.6896;
		                for (int i = 0; i < 4; i++){
		            
		                    //Limiting forces    
		                    //if force is really low just round to 0

		                    if(f[i] < f_min){
		                        f[i] = f_min;
		                    } else if(f[i] > f_max){
		                        f[i] = f_max;
		                    }
							//find exact pwm hardware signal
							pwms[i] = 1100 + ((-b+std::sqrt(b*b - 4*a*(c-f[i])))/(2*a));
		
		                }
		
		                // This function should apply a static mapping from commanded motor
		                // force [N] to PWM [us]. Each PWM channel corresponds to the
		                // appropriately numbered motor, and should be bounded to remain
		                // within 1120 [us] and 1800 [us], and the mapping should be
		                // defined such that a PWM of 1100 [us] would correspond to zero
		                // thrust on any given motor. See the project document for more
		                // details.
			}
	};
}
