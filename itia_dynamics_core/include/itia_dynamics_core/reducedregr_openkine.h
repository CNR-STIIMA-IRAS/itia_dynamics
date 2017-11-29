
 // -------------------------------------------------------------------------------- 
 // Copyright (c) 2017 CNR-ITIA <iras@itia.cnr.it>
 // All rights reserved.
 //
 // Redistribution and use in source and binary forms, with or without
 // modification, are permitted provided that the following conditions are met:
 //
 // 1. Redistributions of source code must retain the above copyright notice,
 // this list of conditions and the following disclaimer.
 // 2. Redistributions in binary form must reproduce the above copyright
 // notice, this list of conditions and the following disclaimer in the
 // documentation and/or other materials provided with the distribution.
 // 3. Neither the name of mosquitto nor the names of its
 // contributors may be used to endorse or promote products derived from
 // this software without specific prior written permission.
 //
 //
 // THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 // AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 // IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 // ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 // LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 // CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 // SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 // INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 // CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 // ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 // POSSIBILITY OF SUCH DAMAGE.
 // -------------------------------------------------------------------------------- 

#ifndef __COMAU_REGRESSOR__
#define __COMAU_REGRESSOR__

# include <Eigen/Core>

namespace itia
{
namespace dynamics
{


inline Eigen::MatrixXd reducedregr_openkine(  const Eigen::VectorXd& qdh,
                                              const Eigen::VectorXd& dqdh,
                                              const Eigen::VectorXd& ddqdh,
                                              const Eigen::VectorXd& g, 
                                              const Eigen::VectorXd& a_DH, 
                                              const Eigen::VectorXd& d_DH, 
                                              const Eigen::MatrixXd& Kr, 
                                              const double& k54, 
                                              const double& k64, 
                                              const double& k65 )
{

  double r01x =  a_DH(0);                                   // L{1}.A; % dati.r01x;  TabellaDH.a(1)
  double r01y =  d_DH(0);                                   //L{1}.D; % dati.r01y;  -TabellaDH.d(1)
  double r12x =  a_DH(1);                                   //L{2}.A; % dati.r12x;  TabellaDH.a(2)
  double r23x =  a_DH(2);                                   //L{3}.A; % dati.r23x;  TabellaDH.a(3)
  double r34y =  d_DH(3);                                   //L{4}.D; % dati.r34y;  TabellaDH.d(4)
  double r56z =  d_DH(5);                                   //L{6}.D; % dati.r56z;  TabellaDH.d(6)

  double gx = g(0);                                         // [m/sec^2]
  double gy = g(1);                                         // [m/sec^2]
  double gz = g(2);                                         // [m/sec^2]



  double q1        = qdh(0); 
  double q2        = qdh(1); 
  double q3        = qdh(2); 
  double q4        = qdh(3); 
  double q5        = qdh(4); 
  double q6        = qdh(5);

  double dq1       = dqdh(0); 
  double dq2       = dqdh(1); 
  double dq3       = dqdh(2); 
  double dq4       = dqdh(3); 
  double dq5       = dqdh(4); 
  double dq6       = dqdh(5);

  double ddq1      = ddqdh(0); 
  double ddq2      = ddqdh(1); 
  double ddq3      = ddqdh(2); 
  double ddq4      = ddqdh(3); 
  double ddq5      = ddqdh(4); 
  double ddq6      = ddqdh(5);



  double w0x   = 0;
  double w0y   = 0;
  double w0z   = 0;

  double w1x   = 0;
  double w1y   = -dq1;
  double w1z   = 0;

  double w2x   = -sin(q2)*dq1;
  double w2y   = -cos(q2)*dq1;
  double w2z   = dq2;

  double w3x   = -dq1*sin(q3+q2);
  double w3y   = -dq2-dq3;
  double w3z   = -dq1*cos(q3+q2);

  double w4x   = -cos(q4)*dq1*sin(q3+q2)+sin(q4)*(-dq2-dq3);
  double w4y   = -dq1*cos(q3+q2)+dq4;
  double w4z   = -sin(q4)*dq1*sin(q3+q2)-cos(q4)*(-dq2-dq3);

  double w5x   = cos(q5)*(-cos(q4)*dq1*sin(q3+q2)+sin(q4)*(-dq2-dq3))+sin(q5)*(-dq1*cos(q3+q2)+dq4);
  double w5y   = sin(q4)*dq1*sin(q3+q2)+cos(q4)*(-dq2-dq3)-dq5;
  double w5z   = -sin(q5)*(-cos(q4)*dq1*sin(q3+q2)+sin(q4)*(-dq2-dq3))+cos(q5)*(-dq1*cos(q3+q2)+dq4);

  double w6x   = cos(q6)*(cos(q5)*(-cos(q4)*dq1*sin(q3+q2)+sin(q4)*(-dq2-dq3))+sin(q5)*(-dq1*cos(q3+q2)+dq4))+sin(q6)*      (sin(q4)*dq1*sin(q3+q2)+cos(q4)*(-dq2-dq3)-dq5);
  double w6y   = -sin(q6)*(cos(q5)*(-cos(q4)*dq1*sin(q3+q2)+sin(q4)*(-dq2-dq3))+sin(q5)*(-dq1*
      cos(q3+q2)+dq4))+cos(q6)*(sin(q4)*dq1*sin(q3+q2)+cos(q4)*(-dq2-dq3)-dq5);
  double w6z   = -sin(q5)*(-cos(q4)*dq1*sin(q3+q2)+sin(q4)*(-dq2-dq3))+cos(q5)*(-dq1*cos(q3+q2)+dq4)+dq6;



  double  dw0x  = 0;
  double dw0y  = 0;
  double dw0z  = 0;

  double dw1x  = 0;
  double dw1y  = -ddq1;
  double dw1z  = 0;

  double dw2x  = cos(q2)*dq2*w1y-sin(q2)*ddq1;
  double dw2y  = -sin(q2)*dq2*w1y-cos(q2)*ddq1;
  double dw2z  = ddq2;

  double dw3x  = dq2*w1y*cos(q3+q2)-ddq1*sin(q3+q2)+cos(q3)*dq3*w2y-sin(q3)*dq3*w2x;
  double dw3y  = -ddq2-ddq3;
  double dw3z  = -dq2*w1y*sin(q3+q2)-ddq1*cos(q3+q2)-sin(q3)*dq3*w2y-cos(q3)*dq3*w2x;

  double dw4x  = cos(q4)*(dq2*w1y*cos(q3+q2)-ddq1*sin(q3+q2)+cos(q3)*dq3*w2y-sin(q3)*dq3*w2x+dq4*w3y)+sin(q4)*(-ddq2-ddq3-dq4*w3x);
  double dw4y  = -dq2*w1y*sin(q3+q2)-ddq1*cos(q3+q2)-sin(q3)*dq3*w2y-cos(q3)*dq3*w2x+ddq4;
  double dw4z  = sin(q4)*(dq2*w1y*cos(q3+q2)-ddq1*sin(q3+q2)+cos(q3)*dq3*w2y-sin(q3)*dq3*w2x+dq4*w3y)-cos(q4)*(-ddq2-ddq3-dq4*w3x);

  double dw5x  = cos(q5)*(cos(q4)*(dq2*w1y*cos(q3+q2)-ddq1*sin(q3+q2)+cos(q3)*dq3*w2y-sin(q3)*dq3*w2x+dq4*      w3y)+sin(q4)*(-ddq2-ddq3-dq4*w3x)+dq5*w4y)+sin(q5)*(-dq2*w1y*sin(q3+q2)-ddq1*cos(q3+q2)-sin(q3)*dq3*
      w2y-cos(q3)*dq3*w2x+ddq4-dq5*w4x);
  double dw5y  = -sin(q4)*(dq2*w1y*cos(q3+q2)-ddq1*sin(q3+q2)+cos(q3)*dq3*w2y-sin(q3)*dq3*w2x+dq4*w3y)+cos(q4)*       (-ddq2-ddq3-dq4*w3x)-ddq5;
  double dw5z  = -sin(q5)*(cos(q4)*(dq2*w1y*cos(q3+q2)-ddq1*sin(q3+q2)+cos(q3)*dq3*w2y-sin(q3)*dq3*w2x+dq4*       w3y)+sin(q4)*(-ddq2-ddq3-dq4*w3x)+dq5*w4y)+cos(q5)*(-dq2*w1y*sin(q3+q2)-ddq1*cos(q3+q2)-sin(q3)*dq3*       w2y-cos(q3)*dq3*w2x+ddq4-dq5*w4x); 

  double dw6x  = cos(q6)*(cos(q5)*(cos(q4)*(dq2*w1y*cos(q3+q2)-ddq1*sin(q3+q2)+cos(q3)*dq3*w2y-sin(q3)*dq3*w2x+dq4*      w3y)+sin(q4)*(-ddq2-ddq3-dq4*w3x)+dq5*w4y)+sin(q5)*(-dq2*w1y*sin(q3+q2)-ddq1*cos(q3+q2)-sin(q3)*dq3*      w2y-cos(q3)*dq3*w2x+ddq4-dq5*w4x)+dq6*w5y)+sin(q6)*(-sin(q4)*(dq2*w1y*cos(q3+q2)-ddq1*      sin(q3+q2)+cos(q3)*dq3*w2y-sin(q3)*dq3*w2x+dq4*w3y)+cos(q4)*(-ddq2-ddq3-dq4*w3x)-ddq5-dq6*w5x);
  double dw6y  = -sin(q6)*(cos(q5)*(cos(q4)*(dq2*w1y*cos(q3+q2)-ddq1*sin(q3+q2)+cos(q3)*dq3*w2y-sin(q3)*dq3*w2x+dq4*      w3y)+sin(q4)*(-ddq2-ddq3-dq4*w3x)+dq5*w4y)+sin(q5)*(-dq2*w1y*sin(q3+q2)-ddq1*cos(q3+q2)-sin(q3)*dq3*      w2y-cos(q3)*dq3*w2x+ddq4-dq5*w4x)+dq6*w5y)+cos(q6)*(-sin(q4)*(dq2*w1y*cos(q3+q2)-ddq1*      sin(q3+q2)+cos(q3)*dq3*w2y-sin(q3)*dq3*w2x+dq4*w3y)+cos(q4)*(-ddq2-ddq3-dq4*w3x)-ddq5-dq6*w5x);
  double dw6z  = -sin(q5)*(cos(q4)*(dq2*w1y*cos(q3+q2)-ddq1*sin(q3+q2)+cos(q3)*dq3*w2y-sin(q3)*dq3*w2x+dq4*      w3y)+sin(q4)*(-ddq2-ddq3-dq4*w3x)+dq5*w4y)+cos(q5)*(-dq2*w1y*sin(q3+q2)-ddq1*cos(q3+q2)-sin(q3)*dq3*      w2y-cos(q3)*dq3*w2x+ddq4-dq5*w4x)+ddq6;


  double  ddp0x = -gx;
  double  ddp0y = -gy;
  double  ddp0z = -gz;

  double ddp1x = -cos(q1)*gx-sin(q1)*gy-pow(w1y,2.0)*r01x;
  double ddp1y = gz;
  double ddp1z = sin(q1)*gx-cos(q1)*gy-dw1y*r01x;

  double ddp2x = cos(q2)*(-cos(q1)*gx-sin(q1)*gy-pow(w1y,2.0)*r01x)+sin(q2)*gz-pow(w2y,2.0)*r12x-pow(w2z,2.0)*r12x;
  double ddp2y = -sin(q2)*(-cos(q1)*gx-sin(q1)*gy-pow(w1y,2.0)*r01x)+cos(q2)*gz+dw2z*r12x+w2x*w2y*r12x;
  double ddp2z = sin(q1)*gx-cos(q1)*gy-dw1y*r01x-dw2y*r12x+w2x*w2z*r12x;

  double ddp3x = cos(q3)*(cos(q2)*(-cos(q1)*gx-sin(q1)*gy-pow(w1y,2.0)*r01x)+sin(q2)*gz-pow(w2y,2.0)*r12x-pow(w2z,2.0)*r12x)+sin(q3)* 
      (-sin(q2)*(-cos(q1)*gx-sin(q1)*gy-pow(w1y,2.0)*r01x)+cos(q2)*gz+dw2z*r12x+w2x*w2y*r12x)-pow(w3y,2.0)*r23x-pow(w3z,2.0)*
      r23x;
  double ddp3y = -sin(q1)*gx+cos(q1)*gy+dw1y*r01x+dw2y*r12x-w2x*w2z*r12x+dw3z*r23x+w3x*w3y*r23x;
  double ddp3z = -sin(q3)*(cos(q2)*(-cos(q1)*gx-sin(q1)*gy-pow(w1y,2.0)*r01x)+sin(q2)*gz-pow(w2y,2.0)*r12x-pow(w2z,2.0)*r12x)+cos(q3)*
      (-sin(q2)*(-cos(q1)*gx-sin(q1)*gy-pow(w1y,2.0)*r01x)+cos(q2)*gz+dw2z*r12x+w2x*w2y*r12x)-dw3y*r23x+w3x*w3z*
      r23x;

  double ddp4x = cos(q4)*(cos(q3)*(cos(q2)*(-cos(q1)*gx-sin(q1)*gy-pow(w1y,2.0)*r01x)+sin(q2)*gz-pow(w2y,2.0)*r12x-pow(w2z,2.0)* 
      r12x)+sin(q3)*(-sin(q2)*(-cos(q1)*gx-sin(q1)*gy-pow(w1y,2.0)*r01x)+cos(q2)*gz+dw2z*r12x+w2x*w2y*r12x)-pow(w3y,2.0)*
      r23x-pow(w3z,2.0)*r23x)+sin(q4)*(-sin(q1)*gx+cos(q1)*gy+dw1y*r01x+dw2y*r12x-w2x*w2z*r12x+dw3z*r23x+w3x*w3y*
      r23x)-dw4z*r34y+w4y*w4x*r34y;
  double ddp4y = -sin(q3)*(cos(q2)*(-cos(q1)*gx-sin(q1)*gy-pow(w1y,2.0)*r01x)+sin(q2)*gz-pow(w2y,2.0)*r12x-pow(w2z,2.0)*r12x)+cos(q3)*
      (-sin(q2)*(-cos(q1)*gx-sin(q1)*gy-pow(w1y,2.0)*r01x)+cos(q2)*gz+dw2z*r12x+w2x*w2y*r12x)-dw3y*r23x+w3x*w3z*
      r23x-pow(w4z,2.0)*r34y-pow(w4x,2.0)*r34y;
  double ddp4z = sin(q4)*(cos(q3)*(cos(q2)*(-cos(q1)*gx-sin(q1)*gy-pow(w1y,2.0)*r01x)+sin(q2)*gz-pow(w2y,2.0)*r12x-pow(w2z,2.0)*
      r12x)+sin(q3)*(-sin(q2)*(-cos(q1)*gx-sin(q1)*gy-pow(w1y,2.0)*r01x)+cos(q2)*gz+dw2z*r12x+w2x*w2y*r12x)-pow(w3y,2.0)*
      r23x-pow(w3z,2.0)*r23x)-cos(q4)*(-sin(q1)*gx+cos(q1)*gy+dw1y*r01x+dw2y*r12x-w2x*w2z*r12x+dw3z*r23x+w3x*w3y*
      r23x)+dw4x*r34y+w4y*w4z*r34y;

  double ddp5x = cos(q5)*(cos(q4)*(cos(q3)*(cos(q2)*(-cos(q1)*gx-sin(q1)*gy-pow(w1y,2.0)*r01x)+sin(q2)*gz-pow(w2y,2.0)*r12x-pow(w2z,2.0)*
      r12x)+sin(q3)*(-sin(q2)*(-cos(q1)*gx-sin(q1)*gy-pow(w1y,2.0)*r01x)+cos(q2)*gz+dw2z*r12x+w2x*w2y*r12x)-pow(w3y,2.0)*
      r23x-pow(w3z,2.0)*r23x)+sin(q4)*(-sin(q1)*gx+cos(q1)*gy+dw1y*r01x+dw2y*r12x-w2x*w2z*r12x+dw3z*r23x+w3x*w3y*
      r23x)-dw4z*r34y+w4y*w4x*r34y)+sin(q5)*(-sin(q3)*(cos(q2)*(-cos(q1)*gx-sin(q1)*gy-pow(w1y,2.0)*r01x)+sin(q2)*
      gz-pow(w2y,2.0)*r12x-pow(w2z,2.0)*r12x)+cos(q3)*(-sin(q2)*(-cos(q1)*gx-sin(q1)*gy-pow(w1y,2.0)*r01x)+cos(q2)*gz+dw2z*
      r12x+w2x*w2y*r12x)-dw3y*r23x+w3x*w3z*r23x-pow(w4z,2.0)*r34y-pow(w4x,2.0)*r34y);
  double ddp5y = -sin(q4)*(cos(q3)*(cos(q2)*(-cos(q1)*gx-sin(q1)*gy-pow(w1y,2.0)*r01x)+sin(q2)*gz-pow(w2y,2.0)*r12x-pow(w2z,2.0)*
      r12x)+sin(q3)*(-sin(q2)*(-cos(q1)*gx-sin(q1)*gy-pow(w1y,2.0)*r01x)+cos(q2)*gz+dw2z*r12x+w2x*w2y*r12x)-pow(w3y,2.0)*
      r23x-pow(w3z,2.0)*r23x)+cos(q4)*(-sin(q1)*gx+cos(q1)*gy+dw1y*r01x+dw2y*r12x-w2x*w2z*r12x+dw3z*r23x+w3x*w3y*
      r23x)-dw4x*r34y-w4y*w4z*r34y;
  double ddp5z = -sin(q5)*(cos(q4)*(cos(q3)*(cos(q2)*(-cos(q1)*gx-sin(q1)*gy-pow(w1y,2.0)*r01x)+sin(q2)*gz-pow(w2y,2.0)*r12x-pow(w2z,2.0)*
      r12x)+sin(q3)*(-sin(q2)*(-cos(q1)*gx-sin(q1)*gy-pow(w1y,2.0)*r01x)+cos(q2)*gz+dw2z*r12x+w2x*w2y*r12x)-pow(w3y,2.0)*
      r23x-pow(w3z,2.0)*r23x)+sin(q4)*(-sin(q1)*gx+cos(q1)*gy+dw1y*r01x+dw2y*r12x-w2x*w2z*r12x+dw3z*r23x+w3x*w3y*
      r23x)-dw4z*r34y+w4y*w4x*r34y)+cos(q5)*(-sin(q3)*(cos(q2)*(-cos(q1)*gx-sin(q1)*gy-pow(w1y,2.0)*r01x)+sin(q2)*
      gz-pow(w2y,2.0)*r12x-pow(w2z,2.0)*r12x)+cos(q3)*(-sin(q2)*(-cos(q1)*gx-sin(q1)*gy-pow(w1y,2.0)*r01x)+cos(q2)*gz+dw2z*
      r12x+w2x*w2y*r12x)-dw3y*r23x+w3x*w3z*r23x-pow(w4z,2.0)*r34y-pow(w4x,2.0)*r34y);

  double ddp6x = cos(q6)*(cos(q5)*(cos(q4)*(cos(q3)*(cos(q2)*(-cos(q1)*gx-sin(q1)*gy-pow(w1y,2.0)*r01x)+sin(q2)*gz-pow(w2y,2.0)*
      r12x-pow(w2z,2.0)*r12x)+sin(q3)*(-sin(q2)*(-cos(q1)*gx-sin(q1)*gy-pow(w1y,2.0)*r01x)+cos(q2)*gz+dw2z*r12x+w2x*w2y*
      r12x)-pow(w3y,2.0)*r23x-pow(w3z,2.0)*r23x)+sin(q4)*(-sin(q1)*gx+cos(q1)*gy+dw1y*r01x+dw2y*r12x-w2x*w2z*r12x+dw3z*
      r23x+w3x*w3y*r23x)-dw4z*r34y+w4y*w4x*r34y)+sin(q5)*(-sin(q3)*(cos(q2)*(-cos(q1)*gx-sin(q1)*gy-pow(w1y,2.0)*
      r01x)+sin(q2)*gz-pow(w2y,2.0)*r12x-pow(w2z,2.0)*r12x)+cos(q3)*(-sin(q2)*(-cos(q1)*gx-sin(q1)*gy-pow(w1y,2.0)*
      r01x)+cos(q2)*gz+dw2z*r12x+w2x*w2y*r12x)-dw3y*r23x+w3x*w3z*r23x-pow(w4z,2.0)*r34y-pow(w4x,2.0)*r34y))+sin(q6)*
      (-sin(q4)*(cos(q3)*(cos(q2)*(-cos(q1)*gx-sin(q1)*gy-pow(w1y,2.0)*r01x)+sin(q2)*gz-pow(w2y,2.0)*r12x-pow(w2z,2.0)*
      r12x)+sin(q3)*(-sin(q2)*(-cos(q1)*gx-sin(q1)*gy-pow(w1y,2.0)*r01x)+cos(q2)*gz+dw2z*r12x+w2x*w2y*r12x)-pow(w3y,2.0)*
      r23x-pow(w3z,2.0)*r23x)+cos(q4)*(-sin(q1)*gx+cos(q1)*gy+dw1y*r01x+dw2y*r12x-w2x*w2z*r12x+dw3z*r23x+w3x*w3y*
      r23x)-dw4x*r34y-w4y*w4z*r34y)+dw6y*r56z+w6z*w6x*r56z;
  double ddp6y = -sin(q6)*(cos(q5)*(cos(q4)*(cos(q3)*(cos(q2)*(-cos(q1)*gx-sin(q1)*gy-pow(w1y,2.0)*r01x)+sin(q2)*gz-pow(w2y,2.0)*
      r12x-pow(w2z,2.0)*r12x)+sin(q3)*(-sin(q2)*(-cos(q1)*gx-sin(q1)*gy-pow(w1y,2.0)*r01x)+cos(q2)*gz+dw2z*r12x+w2x*w2y*
      r12x)-pow(w3y,2.0)*r23x-pow(w3z,2.0)*r23x)+sin(q4)*(-sin(q1)*gx+cos(q1)*gy+dw1y*r01x+dw2y*r12x-w2x*w2z*r12x+dw3z*
      r23x+w3x*w3y*r23x)-dw4z*r34y+w4y*w4x*r34y)+sin(q5)*(-sin(q3)*(cos(q2)*(-cos(q1)*gx-sin(q1)*gy-pow(w1y,2.0)*
      r01x)+sin(q2)*gz-pow(w2y,2.0)*r12x-pow(w2z,2.0)*r12x)+cos(q3)*(-sin(q2)*(-cos(q1)*gx-sin(q1)*gy-pow(w1y,2.0)*
      r01x)+cos(q2)*gz+dw2z*r12x+w2x*w2y*r12x)-dw3y*r23x+w3x*w3z*r23x-pow(w4z,2.0)*r34y-pow(w4x,2.0)*r34y))+cos(q6)*
      (-sin(q4)*(cos(q3)*(cos(q2)*(-cos(q1)*gx-sin(q1)*gy-pow(w1y,2.0)*r01x)+sin(q2)*gz-pow(w2y,2.0)*r12x-pow(w2z,2.0)*
      r12x)+sin(q3)*(-sin(q2)*(-cos(q1)*gx-sin(q1)*gy-pow(w1y,2.0)*r01x)+cos(q2)*gz+dw2z*r12x+w2x*w2y*r12x)-pow(w3y,2.0)*
      r23x-pow(w3z,2.0)*r23x)+cos(q4)*(-sin(q1)*gx+cos(q1)*gy+dw1y*r01x+dw2y*r12x-w2x*w2z*r12x+dw3z*r23x+w3x*w3y*
      r23x)-dw4x*r34y-w4y*w4z*r34y)-dw6x*r56z+w6z*w6y*r56z;
  double ddp6z = -sin(q5)*(cos(q4)*(cos(q3)*(cos(q2)*(-cos(q1)*gx-sin(q1)*gy-pow(w1y,2.0)*r01x)+sin(q2)*gz-pow(w2y,2.0)*r12x-pow(w2z,2.0)*
      r12x)+sin(q3)*(-sin(q2)*(-cos(q1)*gx-sin(q1)*gy-pow(w1y,2.0)*r01x)+cos(q2)*gz+dw2z*r12x+w2x*w2y*r12x)-pow(w3y,2.0)*
      r23x-pow(w3z,2.0)*r23x)+sin(q4)*(-sin(q1)*gx+cos(q1)*gy+dw1y*r01x+dw2y*r12x-w2x*w2z*r12x+dw3z*r23x+w3x*w3y*
      r23x)-dw4z*r34y+w4y*w4x*r34y)+cos(q5)*(-sin(q3)*(cos(q2)*(-cos(q1)*gx-sin(q1)*gy-pow(w1y,2.0)*r01x)+sin(q2)*
      gz-pow(w2y,2.0)*r12x-pow(w2z,2.0)*r12x)+cos(q3)*(-sin(q2)*(-cos(q1)*gx-sin(q1)*gy-pow(w1y,2.0)*r01x)+cos(q2)*gz+dw2z*
      r12x+w2x*w2y*r12x)-dw3y*r23x+w3x*w3z*r23x-pow(w4z,2.0)*r34y-pow(w4x,2.0)*r34y)-pow(w6x,2.0)*r56z-pow(w6y,2.0)*r56z;


  Eigen::Matrix<double, 6, 40> Phi;
  Phi.setZero();
  // Coppia: 1
  Phi(0,0)  = ddp1z*sin(q2)-r01x*(dw2x+w2y*w2z);
  Phi(0,1)  = (dw2y-w2x*w2z)*sin(q2)+(dw2x+w2y*w2z)*cos(q2);
  Phi(0,2)  = (pow(w2y,2.0)-pow(w2z,2.0))*sin(q2)+(dw2z-w2y*w2x)*cos(q2);
  Phi(0,3)  = dw3y*sin(q3+q2)-w3z*w3x*sin(q3+q2)+pow(w3x,2.0)*cos(q3+q2)-pow(w3y,2.0)*cos(q3+q2);
  Phi(0,4)  = pow(w3y,2.0)*sin(q3+q2)-pow(w3z,2.0)*sin(q3+q2)+dw3y*cos(q3+q2)+w3z*w3x*cos(q3+q2);
  Phi(0,5)  = 0;
  Phi(0,6)  = (ddp3z*sin(q4)*cos(q3)-(-sin(q4)*ddp3x+cos(q4)*ddp3y+r23x*(sin(q4)*(-pow(w4y,2.0)-pow(w4z,2.0))-cos(q4)*
              (-dw4y+w4x*w4z)))*sin(q3))*sin(q2)+(ddp3z*sin(q4)*sin(q3)+(-sin(q4)*ddp3x+cos(q4)*ddp3y+r23x*
              (sin(q4)*(-pow(w4y,2.0)-pow(w4z,2.0))-cos(q4)*(-dw4y+w4x*w4z)))*cos(q3)-r12x*(-sin(q4)*
              (-pow(w4y,2.0)-pow(w4z,2.0))+cos(q4)*(-dw4y+w4x*w4z)))*cos(q2)-r01x*(-sin(q4)*(-pow(w4y,2.0)-pow(w4z,2.0))+cos(q4)*
              (-dw4y+w4x*w4z));
  Phi(0,7)  = (((dw4y-w4x*w4z)*cos(q4)+(pow(w4x,2.0)-pow(w4y,2.0))*sin(q4))*cos(q3)-(dw4x+w4y*w4z)*sin(q3))*
              sin(q2)+(((dw4y-w4x*w4z)*cos(q4)+(pow(w4x,2.0)-pow(w4y,2.0))*sin(q4))*sin(q3)+(dw4x+w4y*w4z)*cos(q3))*
              cos(q2);
  Phi(0,8)  = (((dw4z+w4y*w4x)*cos(q4)+(dw4x-w4y*w4z)*sin(q4))*cos(q3)-(pow(w4z,2.0)-pow(w4x,2.0))*sin(q3))*
              sin(q2)+(((dw4z+w4y*w4x)*cos(q4)+(dw4x-w4y*w4z)*sin(q4))*sin(q3)+(pow(w4z,2.0)-pow(w4x,2.0))*cos(q3))*
              cos(q2);
  Phi(0,9) = 0;
  Phi(0,10) = (((ddp4z*sin(q5)+r34y*(-dw5z-w5x*w5y))*cos(q4)+(-sin(q5)*ddp4x+cos(q5)*ddp4y-r34y*(cos(q5)*
              (-pow(w5y,2.0)-pow(w5z,2.0))-sin(q5)*(-dw5y+w5x*w5z)))*sin(q4))*cos(q3)-(-ddp4z*cos(q5)+r23x*(sin(q4)*
              (cos(q5)*(-pow(w5y,2.0)-pow(w5z,2.0))-sin(q5)*(-dw5y+w5x*w5z))-cos(q4)*(-dw5z-w5x*w5y)))*sin(q3))*
              sin(q2)+(((ddp4z*sin(q5)+r34y*(-dw5z-w5x*w5y))*cos(q4)+(-sin(q5)*ddp4x+cos(q5)*ddp4y-r34y*
              (cos(q5)*(-pow(w5y,2.0)-pow(w5z,2.0))-sin(q5)*(-dw5y+w5x*w5z)))*sin(q4))*sin(q3)+(-ddp4z*cos(q5)+r23x*
              (sin(q4)*(cos(q5)*(-pow(w5y,2.0)-pow(w5z,2.0))-sin(q5)*(-dw5y+w5x*w5z))-cos(q4)*(-dw5z-w5x*w5y)))*
              cos(q3)-r12x*(-sin(q4)*(cos(q5)*(-pow(w5y,2.0)-pow(w5z,2.0))-sin(q5)*(-dw5y+w5x*w5z))+cos(q4)*(-dw5z-w5x*
              w5y)))*cos(q2)-r01x*(-sin(q4)*(cos(q5)*(-pow(w5y,2.0)-pow(w5z,2.0))-sin(q5)*(-dw5y+w5x*w5z))+cos(q4)*
              (-dw5z-w5x*w5y));
  Phi(0,11) = ((((dw5y-w5x*w5z)*cos(q5)-(pow(w5x,2.0)-pow(w5y,2.0))*sin(q5))*cos(q4)+(-dw5x-w5y*w5z)*sin(q4))*
              cos(q3)-((dw5y-w5x*w5z)*sin(q5)+(pow(w5x,2.0)-pow(w5y,2.0))*cos(q5))*sin(q3))*sin(q2)+((((dw5y-w5x*w5z)*
              cos(q5)-(pow(w5x,2.0)-pow(w5y,2.0))*sin(q5))*cos(q4)+(-dw5x-w5y*w5z)*sin(q4))*sin(q3)+((dw5y-w5x*w5z)*
              sin(q5)+(pow(w5x,2.0)-pow(w5y,2.0))*cos(q5))*cos(q3))*cos(q2);
  Phi(0,12) = ((((dw5z+w5x*w5y)*cos(q5)-(dw5x-w5y*w5z)*sin(q5))*cos(q4)+(-pow(w5z,2.0)+pow(w5x,2.0))*sin(q4))*
              cos(q3)-((dw5z+w5x*w5y)*sin(q5)+(dw5x-w5y*w5z)*cos(q5))*sin(q3))*sin(q2)+((((dw5z+w5x*w5y)*
              cos(q5)-(dw5x-w5y*w5z)*sin(q5))*cos(q4)+(-pow(w5z,2.0)+pow(w5x,2.0))*sin(q4))*sin(q3)+((dw5z+w5x*w5y)*
              sin(q5)+(dw5x-w5y*w5z)*cos(q5))*cos(q3))*cos(q2);
  Phi(0,13) = ((((pow(w5y,2.0)-pow(w5z,2.0))*cos(q5)-(dw5y+w5x*w5z)*sin(q5))*cos(q4)+(-dw5z+w5x*w5y)*sin(q4))*
              cos(q3)-((pow(w5y,2.0)-pow(w5z,2.0))*sin(q5)+cos(q5)*(dw5y+w5x*w5z))*sin(q3))*sin(q2)+((((pow(w5y,2.0)-pow(w5z,2.0))*
              cos(q5)-(dw5y+w5x*w5z)*sin(q5))*cos(q4)+(-dw5z+w5x*w5y)*sin(q4))*sin(q3)+((pow(w5y,2.0)-pow(w5z,2.0))*
              sin(q5)+cos(q5)*(dw5y+w5x*w5z))*cos(q3))*cos(q2);
  Phi(0,14) = 0;
  Phi(0,15) = (((ddp5z*sin(q6)*cos(q5)-(-sin(q6)*ddp5x+cos(q6)*ddp5y)*sin(q5)+r34y*(-sin(q6)*
              (-pow(w6y,2.0)-pow(w6z,2.0))-cos(q6)*(dw6z+w6y*w6x)))*cos(q4)+(ddp5z*cos(q6)-r34y*(cos(q5)*(cos(q6)*
              (-pow(w6y,2.0)-pow(w6z,2.0))-(dw6z+w6y*w6x)*sin(q6))-sin(q5)*(-dw6y+w6z*w6x)))*sin(q4))*cos(q3)-(ddp5z*
              sin(q6)*sin(q5)+(-sin(q6)*ddp5x+cos(q6)*ddp5y)*cos(q5)+r23x*(sin(q4)*(cos(q5)*(cos(q6)*
              (-pow(w6y,2.0)-pow(w6z,2.0))-(dw6z+w6y*w6x)*sin(q6))-sin(q5)*(-dw6y+w6z*w6x))-cos(q4)*(-sin(q6)*
              (-pow(w6y,2.0)-pow(w6z,2.0))-cos(q6)*(dw6z+w6y*w6x))))*sin(q3))*sin(q2)+(((ddp5z*sin(q6)*cos(q5)-(-sin(q6)*
              ddp5x+cos(q6)*ddp5y)*sin(q5)+r34y*(-sin(q6)*(-pow(w6y,2.0)-pow(w6z,2.0))-cos(q6)*(dw6z+w6y*w6x)))*
              cos(q4)+(ddp5z*cos(q6)-r34y*(cos(q5)*(cos(q6)*(-pow(w6y,2.0)-pow(w6z,2.0))-(dw6z+w6y*w6x)*sin(q6))-sin(q5)*
              (-dw6y+w6z*w6x)))*sin(q4))*sin(q3)+(ddp5z*sin(q6)*sin(q5)+(-sin(q6)*ddp5x+cos(q6)*ddp5y)*
              cos(q5)+r23x*(sin(q4)*(cos(q5)*(cos(q6)*(-pow(w6y,2.0)-pow(w6z,2.0))-(dw6z+w6y*w6x)*sin(q6))-sin(q5)*
              (-dw6y+w6z*w6x))-cos(q4)*(-sin(q6)*(-pow(w6y,2.0)-pow(w6z,2.0))-cos(q6)*(dw6z+w6y*w6x))))*cos(q3)-r12x*
              (-sin(q4)*(cos(q5)*(cos(q6)*(-pow(w6y,2.0)-pow(w6z,2.0))-(dw6z+w6y*w6x)*sin(q6))-sin(q5)*(-dw6y+w6z*
              w6x))+cos(q4)*(-sin(q6)*(-pow(w6y,2.0)-pow(w6z,2.0))-cos(q6)*(dw6z+w6y*w6x))))*cos(q2)-r01x*(-sin(q4)*
              (cos(q5)*(cos(q6)*(-pow(w6y,2.0)-pow(w6z,2.0))-(dw6z+w6y*w6x)*sin(q6))-sin(q5)*(-dw6y+w6z*w6x))+cos(q4)*
              (-sin(q6)*(-pow(w6y,2.0)-pow(w6z,2.0))-cos(q6)*(dw6z+w6y*w6x)));
  Phi(0,16) = (((ddp5z*cos(q6)*cos(q5)-(-cos(q6)*ddp5x-sin(q6)*ddp5y)*sin(q5)+r34y*(-sin(q6)*(-dw6z+w6y*
              w6x)-cos(q6)*(-pow(w6z,2.0)-pow(w6x,2.0))))*cos(q4)+(-ddp5z*sin(q6)-r34y*(cos(q5)*((-dw6z+w6y*w6x)*
              cos(q6)-sin(q6)*(-pow(w6z,2.0)-pow(w6x,2.0)))-sin(q5)*(dw6x+w6z*w6y)))*sin(q4))*cos(q3)-(ddp5z*cos(q6)*
              sin(q5)+(-cos(q6)*ddp5x-sin(q6)*ddp5y)*cos(q5)+r23x*(sin(q4)*(cos(q5)*((-dw6z+w6y*w6x)*
              cos(q6)-sin(q6)*(-pow(w6z,2.0)-pow(w6x,2.0)))-sin(q5)*(dw6x+w6z*w6y))-cos(q4)*(-sin(q6)*(-dw6z+w6y*
              w6x)-cos(q6)*(-pow(w6z,2.0)-pow(w6x,2.0)))))*sin(q3))*sin(q2)+(((ddp5z*cos(q6)*cos(q5)-(-cos(q6)*
              ddp5x-sin(q6)*ddp5y)*sin(q5)+r34y*(-sin(q6)*(-dw6z+w6y*w6x)-cos(q6)*(-pow(w6z,2.0)-pow(w6x,2.0))))*
              cos(q4)+(-ddp5z*sin(q6)-r34y*(cos(q5)*((-dw6z+w6y*w6x)*cos(q6)-sin(q6)*
              (-pow(w6z,2.0)-pow(w6x,2.0)))-sin(q5)*(dw6x+w6z*w6y)))*sin(q4))*sin(q3)+(ddp5z*cos(q6)*sin(q5)+(-cos(q6)*
              ddp5x-sin(q6)*ddp5y)*cos(q5)+r23x*(sin(q4)*(cos(q5)*((-dw6z+w6y*w6x)*cos(q6)-sin(q6)*
              (-pow(w6z,2.0)-pow(w6x,2.0)))-sin(q5)*(dw6x+w6z*w6y))-cos(q4)*(-sin(q6)*(-dw6z+w6y*w6x)-cos(q6)*
              (-pow(w6z,2.0)-pow(w6x,2.0)))))*cos(q3)-r12x*(-sin(q4)*(cos(q5)*((-dw6z+w6y*w6x)*cos(q6)-sin(q6)*
              (-pow(w6z,2.0)-pow(w6x,2.0)))-sin(q5)*(dw6x+w6z*w6y))+cos(q4)*(-sin(q6)*(-dw6z+w6y*w6x)-cos(q6)*
              (-pow(w6z,2.0)-pow(w6x,2.0)))))*cos(q2)-r01x*(-sin(q4)*(cos(q5)*((-dw6z+w6y*w6x)*cos(q6)-sin(q6)*
              (-pow(w6z,2.0)-pow(w6x,2.0)))-sin(q5)*(dw6x+w6z*w6y))+cos(q4)*(-sin(q6)*(-dw6z+w6y*w6x)-cos(q6)*
              (-pow(w6z,2.0)-pow(w6x,2.0))));
  Phi(0,17) = (((((dw6y-w6z*w6x)*cos(q6)-(dw6x+w6z*w6y)*sin(q6))*cos(q5)-(pow(w6x,2.0)-pow(w6y,2.0))*sin(q5))*
              cos(q4)+(-(dw6y-w6z*w6x)*sin(q6)-(dw6x+w6z*w6y)*cos(q6))*sin(q4))*cos(q3)-(((dw6y-w6z*w6x)*
              cos(q6)-(dw6x+w6z*w6y)*sin(q6))*sin(q5)+(pow(w6x,2.0)-pow(w6y,2.0))*cos(q5))*sin(q3))*sin(q2)+(((((dw6y-w6z*
              w6x)*cos(q6)-(dw6x+w6z*w6y)*sin(q6))*cos(q5)-(pow(w6x,2.0)-pow(w6y,2.0))*sin(q5))*cos(q4)+(-(dw6y-w6z*w6x)*
              sin(q6)-(dw6x+w6z*w6y)*cos(q6))*sin(q4))*sin(q3)+(((dw6y-w6z*w6x)*cos(q6)-(dw6x+w6z*w6y)*
              sin(q6))*sin(q5)+(pow(w6x,2.0)-pow(w6y,2.0))*cos(q5))*cos(q3))*cos(q2);
  Phi(0,18) = ((((cos(q6)*(dw6z+w6y*w6x)-(pow(w6z,2.0)-pow(w6x,2.0))*sin(q6))*cos(q5)-(dw6x-w6z*w6y)*sin(q5))*
              cos(q4)+(-(dw6z+w6y*w6x)*sin(q6)-(pow(w6z,2.0)-pow(w6x,2.0))*cos(q6))*sin(q4))*cos(q3)-((cos(q6)*(dw6z+w6y*
              w6x)-(pow(w6z,2.0)-pow(w6x,2.0))*sin(q6))*sin(q5)+(dw6x-w6z*w6y)*cos(q5))*sin(q3))*sin(q2)+((((cos(q6)*
              (dw6z+w6y*w6x)-(pow(w6z,2.0)-pow(w6x,2.0))*sin(q6))*cos(q5)-(dw6x-w6z*w6y)*sin(q5))*cos(q4)+(-(dw6z+w6y*
              w6x)*sin(q6)-(pow(w6z,2.0)-pow(w6x,2.0))*cos(q6))*sin(q4))*sin(q3)+((cos(q6)*(dw6z+w6y*w6x)-(pow(w6z,2.0)-pow(w6x,2.0))*
              sin(q6))*sin(q5)+(dw6x-w6z*w6y)*cos(q5))*cos(q3))*cos(q2);
  Phi(0,19) = (((((pow(w6y,2.0)-pow(w6z,2.0))*cos(q6)-(dw6z-w6y*w6x)*sin(q6))*cos(q5)-(dw6y+w6z*w6x)*sin(q5))*
              cos(q4)+(-(pow(w6y,2.0)-pow(w6z,2.0))*sin(q6)-(dw6z-w6y*w6x)*cos(q6))*sin(q4))*cos(q3)-(((pow(w6y,2.0)-pow(w6z,2.0))*
              cos(q6)-(dw6z-w6y*w6x)*sin(q6))*sin(q5)+(dw6y+w6z*w6x)*cos(q5))*sin(q3))*
              sin(q2)+(((((pow(w6y,2.0)-pow(w6z,2.0))*cos(q6)-(dw6z-w6y*w6x)*sin(q6))*cos(q5)-(dw6y+w6z*w6x)*sin(q5))*
              cos(q4)+(-(pow(w6y,2.0)-pow(w6z,2.0))*sin(q6)-(dw6z-w6y*w6x)*cos(q6))*sin(q4))*sin(q3)+(((pow(w6y,2.0)-pow(w6z,2.0))*
              cos(q6)-(dw6z-w6y*w6x)*sin(q6))*sin(q5)+(dw6y+w6z*w6x)*cos(q5))*cos(q3))*cos(q2);
  Phi(0,20) = ((((cos(q6)*w6z*w6y+sin(q6)*w6z*w6x)*cos(q5)-dw6z*sin(q5))*cos(q4)+(-w6z*w6y*sin(q6)+w6z*w6x*
              cos(q6))*sin(q4))*cos(q3)-((cos(q6)*w6z*w6y+sin(q6)*w6z*w6x)*sin(q5)+dw6z*cos(q5))*sin(q3))*
              sin(q2)+((((cos(q6)*w6z*w6y+sin(q6)*w6z*w6x)*cos(q5)-dw6z*sin(q5))*cos(q4)+(-w6z*w6y*
              sin(q6)+w6z*w6x*cos(q6))*sin(q4))*sin(q3)+((cos(q6)*w6z*w6y+sin(q6)*w6z*w6x)*sin(q5)+dw6z*
              cos(q5))*cos(q3))*cos(q2);
  Phi(0,21) = 0;
  Phi(0,22) = dw1y;
  Phi(0,23) = -ddp1z*cos(q2)-r01x*(-dw2y+w2x*w2z);
  Phi(0,24) = dw2x*sin(q2)+w2x*w2z*cos(q2);
  Phi(0,25) = (dw2z+w2y*w2x)*sin(q2)+(pow(w2z,2.0)-pow(w2x,2.0))*cos(q2);
  Phi(0,26) = w2z*(w2y*sin(q2)-w2x*cos(q2));
  Phi(0,27) = -ddp2z*cos(q3+q2)+cos(q2)*r12x*dw3z+cos(q2)*r12x*w3x*w3y+r01x*dw3z+r01x*w3x*w3y;
  Phi(0,28) = ddp2z*sin(q3+q2)-cos(q2)*r12x*dw3x+cos(q2)*r12x*w3z*w3y-r01x*dw3x+r01x*w3z*w3y;
  Phi(0,29) = dw3x*sin(q3+q2)-w3x*w3y*cos(q3+q2);
  Phi(0,30) = dw3z*sin(q3+q2)+w3y*w3x*sin(q3+q2)+dw3x*cos(q3+q2)-w3y*w3z*cos(q3+q2);
  Phi(0,31) = -w3y*w3z*sin(q3+q2)+w3x*w3y*cos(q3+q2);
  Phi(0,32) = (-ddp3z*cos(q4)*cos(q3)-(cos(q4)*ddp3x+sin(q4)*ddp3y+r23x*(sin(q4)*(dw4y+w4x*w4z)-cos(q4)*
              (-pow(w4x,2.0)-pow(w4y,2.0))))*sin(q3))*sin(q2)+(-ddp3z*cos(q4)*sin(q3)+(cos(q4)*ddp3x+sin(q4)*ddp3y+r23x*
              (sin(q4)*(dw4y+w4x*w4z)-cos(q4)*(-pow(w4x,2.0)-pow(w4y,2.0))))*cos(q3)-r12x*(-sin(q4)*(dw4y+w4x*
              w4z)+cos(q4)*(-pow(w4x,2.0)-pow(w4y,2.0))))*cos(q2)-r01x*(-sin(q4)*(dw4y+w4x*w4z)+cos(q4)*(-pow(w4x,2.0)-pow(w4y,2.0)));
  Phi(0,33) = 0.5*dw4x*sin(q2-q4+q3)+0.5*dw4x*sin(q2+q4+q3)+0.5*w4y*w4x*cos(q2+q4+q3)-0.5*w4y*w4x*
              cos(q2-q4+q3)+w4z*w4x*cos(q3+q2);
  Phi(0,34) = -0.5*w4z*w4y*sin(q2-q4+q3)-0.5*w4z*w4y*sin(q2+q4+q3)-0.5*w4y*w4x*cos(q2+q4+q3)+0.5*w4y*w4x*
              cos(q2-q4+q3)+dw4y*cos(q3+q2);
  Phi(0,35) = (((pow(w4y,2.0)-pow(w4z,2.0))*cos(q4)+sin(q4)*(dw4y+w4x*w4z))*cos(q3)-(dw4z-w4y*w4x)*sin(q3))*
              sin(q2)+(((pow(w4y,2.0)-pow(w4z,2.0))*cos(q4)+sin(q4)*(dw4y+w4x*w4z))*sin(q3)+(dw4z-w4y*w4x)*cos(q3))*
              cos(q2);
  Phi(0,36) = (((ddp4z*cos(q5)+r34y*(dw5x-w5y*w5z))*cos(q4)+(-cos(q5)*ddp4x-sin(q5)*ddp4y-r34y*(cos(q5)*
              (dw5y+w5x*w5z)-sin(q5)*(-pow(w5x,2.0)-pow(w5y,2.0))))*sin(q4))*cos(q3)-(ddp4z*sin(q5)+r23x*(sin(q4)*
              (cos(q5)*(dw5y+w5x*w5z)-sin(q5)*(-pow(w5x,2.0)-pow(w5y,2.0)))-cos(q4)*(dw5x-w5y*w5z)))*sin(q3))*
              sin(q2)+(((ddp4z*cos(q5)+r34y*(dw5x-w5y*w5z))*cos(q4)+(-cos(q5)*ddp4x-sin(q5)*ddp4y-r34y*
              (cos(q5)*(dw5y+w5x*w5z)-sin(q5)*(-pow(w5x,2.0)-pow(w5y,2.0))))*sin(q4))*sin(q3)+(ddp4z*sin(q5)+r23x*
              (sin(q4)*(cos(q5)*(dw5y+w5x*w5z)-sin(q5)*(-pow(w5x,2.0)-pow(w5y,2.0)))-cos(q4)*(dw5x-w5y*w5z)))*
              cos(q3)-r12x*(-sin(q4)*(cos(q5)*(dw5y+w5x*w5z)-sin(q5)*(-pow(w5x,2.0)-pow(w5y,2.0)))+cos(q4)*(dw5x-w5y*
              w5z)))*cos(q2)-r01x*(-sin(q4)*(cos(q5)*(dw5y+w5x*w5z)-sin(q5)*(-pow(w5x,2.0)-pow(w5y,2.0)))+cos(q4)*
              (dw5x-w5y*w5z));
  Phi(0,37) = (((dw5x*cos(q5)+w5y*w5x*sin(q5))*cos(q4)-w5z*w5x*sin(q4))*cos(q3)-(dw5x*sin(q5)-w5y*w5x*
              cos(q5))*sin(q3))*sin(q2)+(((dw5x*cos(q5)+w5y*w5x*sin(q5))*cos(q4)-w5z*w5x*sin(q4))*
              sin(q3)+(dw5x*sin(q5)-w5y*w5x*cos(q5))*cos(q3))*cos(q2);
  Phi(0,38) = (((-w5y*w5z*cos(q5)-w5y*w5x*sin(q5))*cos(q4)-dw5y*sin(q4))*cos(q3)-(-w5y*w5z*sin(q5)+w5y*w5x*
              cos(q5))*sin(q3))*sin(q2)+(((-w5y*w5z*cos(q5)-w5y*w5x*sin(q5))*cos(q4)-dw5y*sin(q4))*
              sin(q3)+(-w5y*w5z*sin(q5)+w5y*w5x*cos(q5))*cos(q3))*cos(q2);
  Phi(0,39) = ((((cos(q6)*dw6x-sin(q6)*w6z*w6x)*cos(q5)+w6y*w6x*sin(q5))*cos(q4)+(-dw6x*sin(q6)-w6z*w6x*
              cos(q6))*sin(q4))*cos(q3)-((cos(q6)*dw6x-sin(q6)*w6z*w6x)*sin(q5)-w6y*w6x*cos(q5))*sin(q3))*
              sin(q2)+((((cos(q6)*dw6x-sin(q6)*w6z*w6x)*cos(q5)+w6y*w6x*sin(q5))*cos(q4)+(-dw6x*sin(q6)-w6z*
              w6x*cos(q6))*sin(q4))*sin(q3)+((cos(q6)*dw6x-sin(q6)*w6z*w6x)*sin(q5)-w6y*w6x*cos(q5))*
              cos(q3))*cos(q2);
  // Coppia: 2
  Phi(1,0)  = -cos(q2)*ddp1x-sin(q2)*ddp1y;
  Phi(1,1)  = pow(w2x,2.0)-pow(w2y,2.0);
  Phi(1,2)  = dw2y+w2x*w2z;
  Phi(1,3)  = -dw3x-w3y*w3z;
  Phi(1,4)  = -dw3z+w3y*w3x;
  Phi(1,5)  = 0;
  Phi(1,6)  = ddp3z*cos(q4)+r23x*(dw4z+w4x*w4y)+r12x*(sin(q3)*(cos(q4)*(-pow(w4y,2.0)-pow(w4z,2.0))+sin(q4)*(-dw4y+w4x*
              w4z))+cos(q3)*(dw4z+w4x*w4y));
  Phi(1,7)  = -sin(q4)*(dw4y-w4x*w4z)+(pow(w4x,2.0)-pow(w4y,2.0))*cos(q4);
  Phi(1,8)  = -(dw4z+w4y*w4x)*sin(q4)+(dw4x-w4y*w4z)*cos(q4);
  Phi(1,9) = 0;
  Phi(1,10) = -(ddp4z*sin(q5)+r34y*(-dw5z-w5x*w5y))*sin(q4)+(-sin(q5)*ddp4x+cos(q5)*ddp4y-r34y*(cos(q5)*
              (-pow(w5y,2.0)-pow(w5z,2.0))-sin(q5)*(-dw5y+w5x*w5z)))*cos(q4)+r23x*(sin(q5)*(-pow(w5y,2.0)-pow(w5z,2.0))+cos(q5)*
              (-dw5y+w5x*w5z))+r12x*(sin(q3)*(cos(q4)*(cos(q5)*(-pow(w5y,2.0)-pow(w5z,2.0))-sin(q5)*(-dw5y+w5x*
              w5z))+sin(q4)*(-dw5z-w5x*w5y))+cos(q3)*(sin(q5)*(-pow(w5y,2.0)-pow(w5z,2.0))+cos(q5)*(-dw5y+w5x*w5z)));
  Phi(1,11) = -((dw5y-w5z*w5x)*cos(q5)-(pow(w5x,2.0)-pow(w5y,2.0))*sin(q5))*sin(q4)+(-dw5x-w5y*w5z)*cos(q4);
  Phi(1,12) = -((dw5z+w5y*w5x)*cos(q5)-(dw5x-w5y*w5z)*sin(q5))*sin(q4)+(-pow(w5z,2.0)+pow(w5x,2.0))*cos(q4);
  Phi(1,13) = -((pow(w5y,2.0)-pow(w5z,2.0))*cos(q5)-(dw5y+w5z*w5x)*sin(q5))*sin(q4)+(-dw5z+w5y*w5x)*cos(q4);
  Phi(1,14) = 0;
  Phi(1,15) = -(ddp5z*sin(q6)*cos(q5)-(-sin(q6)*ddp5x+cos(q6)*ddp5y)*sin(q5)+r34y*(-sin(q6)*
              (-pow(w6y,2.0)-pow(w6z,2.0))-cos(q6)*(dw6z+w6x*w6y)))*sin(q4)+(ddp5z*cos(q6)-r34y*(cos(q5)*(cos(q6)*
              (-pow(w6y,2.0)-pow(w6z,2.0))-sin(q6)*(dw6z+w6x*w6y))-sin(q5)*(-dw6y+w6x*w6z)))*cos(q4)+r23x*(sin(q5)*
              (cos(q6)*(-pow(w6y,2.0)-pow(w6z,2.0))-sin(q6)*(dw6z+w6x*w6y))+cos(q5)*(-dw6y+w6x*w6z))+r12x*(sin(q3)*
              (cos(q4)*(cos(q5)*(cos(q6)*(-pow(w6y,2.0)-pow(w6z,2.0))-sin(q6)*(dw6z+w6x*w6y))-sin(q5)*(-dw6y+w6x*
              w6z))+sin(q4)*(-sin(q6)*(-pow(w6y,2.0)-pow(w6z,2.0))-cos(q6)*(dw6z+w6x*w6y)))+cos(q3)*(sin(q5)*(cos(q6)*
              (-pow(w6y,2.0)-pow(w6z,2.0))-sin(q6)*(dw6z+w6x*w6y))+cos(q5)*(-dw6y+w6x*w6z)));
  Phi(1,16) = -(ddp5z*cos(q6)*cos(q5)-(-cos(q6)*ddp5x-sin(q6)*ddp5y)*sin(q5)+r34y*(-sin(q6)*(-dw6z+w6x*
              w6y)-cos(q6)*(-pow(w6z,2.0)-pow(w6x,2.0))))*sin(q4)+(-ddp5z*sin(q6)-r34y*(cos(q5)*(cos(q6)*(-dw6z+w6x*
              w6y)-sin(q6)*(-pow(w6z,2.0)-pow(w6x,2.0)))-sin(q5)*(dw6x+w6y*w6z)))*cos(q4)+r23x*(sin(q5)*(cos(q6)*
              (-dw6z+w6x*w6y)-sin(q6)*(-pow(w6z,2.0)-pow(w6x,2.0)))+cos(q5)*(dw6x+w6y*w6z))+r12x*(sin(q3)*(cos(q4)*
              (cos(q5)*(cos(q6)*(-dw6z+w6x*w6y)-sin(q6)*(-pow(w6z,2.0)-pow(w6x,2.0)))-sin(q5)*(dw6x+w6y*w6z))+sin(q4)*
              (-sin(q6)*(-dw6z+w6x*w6y)-cos(q6)*(-pow(w6z,2.0)-pow(w6x,2.0))))+cos(q3)*(sin(q5)*(cos(q6)*(-dw6z+w6x*
              w6y)-sin(q6)*(-pow(w6z,2.0)-pow(w6x,2.0)))+cos(q5)*(dw6x+w6y*w6z)));
  Phi(1,17) = -(((dw6y-w6z*w6x)*cos(q6)-(dw6x+w6y*w6z)*sin(q6))*cos(q5)-(pow(w6x,2.0)-pow(w6y,2.0))*sin(q5))*
              sin(q4)+(-(dw6y-w6z*w6x)*sin(q6)-(dw6x+w6y*w6z)*cos(q6))*cos(q4);
  Phi(1,18) = -(((dw6z+w6y*w6x)*cos(q6)-(pow(w6z,2.0)-pow(w6x,2.0))*sin(q6))*cos(q5)-(dw6x-w6y*w6z)*sin(q5))*
              sin(q4)+(-(dw6z+w6y*w6x)*sin(q6)-(pow(w6z,2.0)-pow(w6x,2.0))*cos(q6))*cos(q4);
  Phi(1,19) = -(((pow(w6y,2.0)-pow(w6z,2.0))*cos(q6)-(dw6z-w6y*w6x)*sin(q6))*cos(q5)-(dw6y+w6x*w6z)*sin(q5))*
              sin(q4)+(-(pow(w6y,2.0)-pow(w6z,2.0))*sin(q6)-(dw6z-w6y*w6x)*cos(q6))*cos(q4);
  Phi(1,20) = -((cos(q6)*w6y*w6z+w6x*w6z*sin(q6))*cos(q5)-dw6z*sin(q5))*sin(q4)+(-sin(q6)*w6y*w6z+w6x*w6z*
              cos(q6))*cos(q4);
  Phi(1,21) = 0;
  Phi(1,22) = 0;
  Phi(1,23) = -sin(q2)*ddp1x+cos(q2)*ddp1y;
  Phi(1,24) = -w2y*w2x;
  Phi(1,25) = dw2x-w2y*w2z;
  Phi(1,26) = dw2z;
  Phi(1,27) = -sin(q3)*ddp2x+cos(q3)*ddp2y+r12x*(sin(q3)*(-pow(w3y,2.0)-pow(w3z,2.0))+cos(q3)*(-dw3y+w3x*w3z));
  Phi(1,28) = -cos(q3)*ddp2x-sin(q3)*ddp2y+r12x*(sin(q3)*(dw3y+w3x*w3z)+cos(q3)*(-pow(w3x,2.0)-pow(w3y,2.0)));
  Phi(1,29) = -w3x*w3z;
  Phi(1,30) = -pow(w3z,2.0)+pow(w3x,2.0);
  Phi(1,31) = -dw3y;
  Phi(1,32) = ddp3z*sin(q4)+r23x*(-dw4x+w4y*w4z)+r12x*(sin(q3)*(cos(q4)*(dw4y+w4x*w4z)+sin(q4)*
              (-pow(w4x,2.0)-pow(w4y,2.0)))+cos(q3)*(-dw4x+w4y*w4z));
  Phi(1,33) = -sin(q4)*dw4x-cos(q4)*w4x*w4y;
  Phi(1,34) = w4y*(sin(q4)*w4z+w4x*cos(q4));
  Phi(1,35) = -(pow(w4y,2.0)-pow(w4z,2.0))*sin(q4)+(dw4y+w4x*w4z)*cos(q4);
  Phi(1,36) = -(ddp4z*cos(q5)+r34y*(dw5x-w5y*w5z))*sin(q4)+(-cos(q5)*ddp4x-sin(q5)*ddp4y-r34y*(cos(q5)*
              (dw5y+w5z*w5x)-sin(q5)*(-pow(w5x,2.0)-pow(w5y,2.0))))*cos(q4)+r23x*(sin(q5)*(dw5y+w5z*w5x)+cos(q5)*
              (-pow(w5x,2.0)-pow(w5y,2.0)))+r12x*(sin(q3)*(cos(q4)*(cos(q5)*(dw5y+w5z*w5x)-sin(q5)*
              (-pow(w5x,2.0)-pow(w5y,2.0)))+sin(q4)*(dw5x-w5y*w5z))+cos(q3)*(sin(q5)*(dw5y+w5z*w5x)+cos(q5)*
              (-pow(w5x,2.0)-pow(w5y,2.0))));
  Phi(1,37) = -(dw5x*cos(q5)+w5y*w5x*sin(q5))*sin(q4)-w5z*w5x*cos(q4);
  Phi(1,38) = -(-w5y*w5z*cos(q5)-w5y*w5x*sin(q5))*sin(q4)-dw5y*cos(q4);
  Phi(1,39) = -((dw6x*cos(q6)-w6z*w6x*sin(q6))*cos(q5)+w6y*w6x*sin(q5))*sin(q4)+(-dw6x*sin(q6)-w6z*w6x*
              cos(q6))*cos(q4);
  // Coppia: 3
  Phi(2,0)  = 0;
  Phi(2,1)  = 0;
  Phi(2,2)  = 0;
  Phi(2,3)  = dw3x+w3z*w3y;
  Phi(2,4)  = dw3z-w3x*w3y;
  Phi(2,5)  = -ddq3;
  Phi(2,6)  = -ddp3z*cos(q4)-r23x*(dw4z+w4y*w4x);
  Phi(2,7)  = (dw4y-w4x*w4z)*sin(q4)-(pow(w4x,2.0)-pow(w4y,2.0))*cos(q4);
  Phi(2,8)  = (dw4z+w4y*w4x)*sin(q4)-(dw4x-w4y*w4z)*cos(q4);
  Phi(2,9) = 0;
  Phi(2,10) = (ddp4z*sin(q5)+r34y*(-dw5z-w5x*w5y))*sin(q4)-(-sin(q5)*ddp4x+cos(q5)*ddp4y-r34y*(cos(q5)*
              (-pow(w5y,2.0)-pow(w5z,2.0))-sin(q5)*(-dw5y+w5x*w5z)))*cos(q4)-r23x*(sin(q5)*(-pow(w5y,2.0)-pow(w5z,2.0))+cos(q5)*
              (-dw5y+w5x*w5z));
  Phi(2,11) = ((dw5y-w5x*w5z)*cos(q5)-(pow(w5x,2.0)-pow(w5y,2.0))*sin(q5))*sin(q4)-(-dw5x-w5y*w5z)*cos(q4);
  Phi(2,12) = ((dw5z+w5x*w5y)*cos(q5)-(dw5x-w5y*w5z)*sin(q5))*sin(q4)-(-pow(w5z,2.0)+pow(w5x,2.0))*cos(q4);
  Phi(2,13) = ((pow(w5y,2.0)-pow(w5z,2.0))*cos(q5)-(dw5y+w5x*w5z)*sin(q5))*sin(q4)-(-dw5z+w5x*w5y)*cos(q4);
  Phi(2,14) = 0;
  Phi(2,15) = (ddp5z*sin(q6)*cos(q5)-(-sin(q6)*ddp5x+cos(q6)*ddp5y)*sin(q5)+r34y*(-sin(q6)*
              (-pow(w6y,2.0)-pow(w6z,2.0))-cos(q6)*(dw6z+w6y*w6x)))*sin(q4)-(ddp5z*cos(q6)-r34y*(cos(q5)*(cos(q6)*
              (-pow(w6y,2.0)-pow(w6z,2.0))-(dw6z+w6y*w6x)*sin(q6))-sin(q5)*(-dw6y+w6z*w6x)))*cos(q4)-r23x*(sin(q5)*
              (cos(q6)*(-pow(w6y,2.0)-pow(w6z,2.0))-(dw6z+w6y*w6x)*sin(q6))+cos(q5)*(-dw6y+w6z*w6x));
  Phi(2,16) = (ddp5z*cos(q6)*cos(q5)-(-cos(q6)*ddp5x-sin(q6)*ddp5y)*sin(q5)+r34y*(-sin(q6)*(-dw6z+w6y*
              w6x)-cos(q6)*(-pow(w6z,2.0)-pow(w6x,2.0))))*sin(q4)-(-ddp5z*sin(q6)-r34y*(cos(q5)*((-dw6z+w6y*w6x)*
              cos(q6)-sin(q6)*(-pow(w6z,2.0)-pow(w6x,2.0)))-sin(q5)*(dw6x+w6z*w6y)))*cos(q4)-r23x*(sin(q5)*((-dw6z+w6y*
              w6x)*cos(q6)-sin(q6)*(-pow(w6z,2.0)-pow(w6x,2.0)))+cos(q5)*(dw6x+w6z*w6y));
  Phi(2,17) = (((dw6y-w6z*w6x)*cos(q6)-(dw6x+w6z*w6y)*sin(q6))*cos(q5)-(pow(w6x,2.0)-pow(w6y,2.0))*sin(q5))*
              sin(q4)-(-(dw6y-w6z*w6x)*sin(q6)-(dw6x+w6z*w6y)*cos(q6))*cos(q4);
  Phi(2,18) = ((cos(q6)*(dw6z+w6y*w6x)-(pow(w6z,2.0)-pow(w6x,2.0))*sin(q6))*cos(q5)-(dw6x-w6z*w6y)*sin(q5))*
              sin(q4)-(-(dw6z+w6y*w6x)*sin(q6)-(pow(w6z,2.0)-pow(w6x,2.0))*cos(q6))*cos(q4);
  Phi(2,19) = (((pow(w6y,2.0)-pow(w6z,2.0))*cos(q6)-(dw6z-w6y*w6x)*sin(q6))*cos(q5)-(dw6y+w6z*w6x)*sin(q5))*
              sin(q4)-(-(pow(w6y,2.0)-pow(w6z,2.0))*sin(q6)-(dw6z-w6y*w6x)*cos(q6))*cos(q4);
  Phi(2,20) = ((cos(q6)*w6z*w6y+sin(q6)*w6z*w6x)*cos(q5)-dw6z*sin(q5))*sin(q4)-(-w6z*w6y*sin(q6)+w6z*w6x*
              cos(q6))*cos(q4);
  Phi(2,21) = 0;
  Phi(2,22) = 0;
  Phi(2,23) = 0;
  Phi(2,24) = 0;
  Phi(2,25) = 0;
  Phi(2,26) = 0;
  Phi(2,27) = sin(q3)*ddp2x-cos(q3)*ddp2y;
  Phi(2,28) = cos(q3)*ddp2x+sin(q3)*ddp2y;
  Phi(2,29) = w3z*w3x;
  Phi(2,30) = pow(w3z,2.0)-pow(w3x,2.0);
  Phi(2,31) = dw3y;
  Phi(2,32) = -ddp3z*sin(q4)-r23x*(-dw4x+w4y*w4z);
  Phi(2,33) = dw4x*sin(q4)+w4x*w4y*cos(q4);
  Phi(2,34) = -w4y*(sin(q4)*w4z+w4x*cos(q4));
  Phi(2,35) = (pow(w4y,2.0)-pow(w4z,2.0))*sin(q4)-(dw4y+w4x*w4z)*cos(q4);
  Phi(2,36) = (ddp4z*cos(q5)+r34y*(dw5x-w5y*w5z))*sin(q4)-(-cos(q5)*ddp4x-sin(q5)*ddp4y-r34y*(cos(q5)*
              (dw5y+w5x*w5z)-sin(q5)*(-pow(w5x,2.0)-pow(w5y,2.0))))*cos(q4)-r23x*((dw5y+w5x*w5z)*sin(q5)+cos(q5)*
              (-pow(w5x,2.0)-pow(w5y,2.0)));
  Phi(2,37) = (dw5x*cos(q5)+w5y*w5x*sin(q5))*sin(q4)+w5x*w5z*cos(q4);
  Phi(2,38) = (-w5y*w5z*cos(q5)-w5y*w5x*sin(q5))*sin(q4)+dw5y*cos(q4);
  Phi(2,39) = ((cos(q6)*dw6x-sin(q6)*w6z*w6x)*cos(q5)+w6y*w6x*sin(q5))*sin(q4)-(-dw6x*sin(q6)-w6z*w6x*
              cos(q6))*cos(q4);
  // Coppia: 4
  Phi(3,0)  = 0;
  Phi(3,1)  = 0;
  Phi(3,2)  = 0;
  Phi(3,3)  = 0;
  Phi(3,4)  = 0;
  Phi(3,5)  = 0;
  Phi(3,6)  = sin(q4)*ddp3x-cos(q4)*ddp3y;
  Phi(3,7)  = -dw4x-w4y*w4z;
  Phi(3,8)  = -pow(w4z,2.0)+pow(w4x,2.0);
  Phi(3,9) = -ddq4;
  Phi(3,10) = ddp4z*cos(q5)+k54*(-sin(q5)*ddp4x+cos(q5)*ddp4y);
  Phi(3,11) = -(dw5y-w5x*w5z)*sin(q5)-(pow(w5x,2.0)-pow(w5y,2.0))*cos(q5)+k54*(-dw5x-w5y*w5z);
  Phi(3,12) = -(dw5z+w5x*w5y)*sin(q5)-(dw5x-w5y*w5z)*cos(q5)+k54*(-pow(w5z,2.0)+pow(w5x,2.0));
  Phi(3,13) = -(pow(w5y,2.0)-pow(w5z,2.0))*sin(q5)-cos(q5)*(dw5y+w5x*w5z)+k54*(-dw5z+w5x*w5y);
  Phi(3,14) = 0;
  Phi(3,15) = -ddp5z*sin(q6)*sin(q5)-(-sin(q6)*ddp5x+cos(q6)*ddp5y)*cos(q5)+k54*ddp5z*cos(q6)-k64*(-sin(q6)*
              ddp5x+cos(q6)*ddp5y);
  Phi(3,16) = -ddp5z*cos(q6)*sin(q5)-(-cos(q6)*ddp5x-sin(q6)*ddp5y)*cos(q5)-k54*ddp5z*sin(q6)-k64*(-cos(q6)*
              ddp5x-sin(q6)*ddp5y);
  Phi(3,17) = -((dw6y-w6z*w6x)*cos(q6)-(dw6x+w6z*w6y)*sin(q6))*sin(q5)-(pow(w6x,2.0)-pow(w6y,2.0))*cos(q5)+k54*
              (-(dw6y-w6z*w6x)*sin(q6)-(dw6x+w6z*w6y)*cos(q6))-k64*(pow(w6x,2.0)-pow(w6y,2.0));
  Phi(3,18) = -(cos(q6)*(dw6z+w6y*w6x)-(pow(w6z,2.0)-pow(w6x,2.0))*sin(q6))*sin(q5)-(dw6x-w6z*w6y)*cos(q5)+k54*
              (-(dw6z+w6y*w6x)*sin(q6)-(pow(w6z,2.0)-pow(w6x,2.0))*cos(q6))-k64*(dw6x-w6z*w6y);
  Phi(3,19) = -((pow(w6y,2.0)-pow(w6z,2.0))*cos(q6)-(dw6z-w6y*w6x)*sin(q6))*sin(q5)-(dw6y+w6z*w6x)*cos(q5)+k54*
              (-(pow(w6y,2.0)-pow(w6z,2.0))*sin(q6)-(dw6z-w6y*w6x)*cos(q6))-k64*(dw6y+w6z*w6x);
  Phi(3,20) = -(cos(q6)*w6z*w6y+sin(q6)*w6z*w6x)*sin(q5)-dw6z*cos(q5)+k54*w6z*(-w6y*sin(q6)+w6x*
              cos(q6))-k64*dw6z;
  Phi(3,21) = 0;
  Phi(3,22) = 0;
  Phi(3,23) = 0;
  Phi(3,24) = 0;
  Phi(3,25) = 0;
  Phi(3,26) = 0;
  Phi(3,27) = 0;
  Phi(3,28) = 0;
  Phi(3,29) = 0;
  Phi(3,30) = 0;
  Phi(3,31) = 0;
  Phi(3,32) = -cos(q4)*ddp3x-sin(q4)*ddp3y;
  Phi(3,33) = -w4x*w4z;
  Phi(3,34) = -dw4y;
  Phi(3,35) = -dw4z+w4y*w4x;
  Phi(3,36) = -ddp4z*sin(q5)+k54*(-cos(q5)*ddp4x-sin(q5)*ddp4y);
  Phi(3,37) = -dw5x*sin(q5)+w5y*w5x*cos(q5)-k54*w5x*w5z;
  Phi(3,38) = -w5y*(-w5z*sin(q5)+w5x*cos(q5))-k54*dw5y;
  Phi(3,39) = -(cos(q6)*dw6x-sin(q6)*w6z*w6x)*sin(q5)+w6y*w6x*cos(q5)+k54*(-dw6x*sin(q6)-w6z*w6x*
              cos(q6))+k64*w6y*w6x;
  // Coppia: 5
  Phi(4,0)  = 0;
  Phi(4,1)  = 0;
  Phi(4,2)  = 0;
  Phi(4,3)  = 0;
  Phi(4,4)  = 0;
  Phi(4,5)  = 0;
  Phi(4,6)  = 0;
  Phi(4,7)  = 0;
  Phi(4,8)  = 0;
  Phi(4,9) = 0;
  Phi(4,10) = -sin(q5)*ddp4x+cos(q5)*ddp4y;
  Phi(4,11) = -dw5x-w5y*w5z;
  Phi(4,12) = -pow(w5z,2.0)+pow(w5x,2.0);
  Phi(4,13) = -dw5z+w5x*w5y;
  Phi(4,14) = k54*ddq4+ddq5;
  Phi(4,15) = ddp5z*cos(q6)-k65*(-sin(q6)*ddp5x+cos(q6)*ddp5y);
  Phi(4,16) = -ddp5z*sin(q6)-k65*(-cos(q6)*ddp5x-sin(q6)*ddp5y);
  Phi(4,17) = -(dw6y-w6z*w6x)*sin(q6)-(dw6x+w6z*w6y)*cos(q6)-k65*(pow(w6x,2.0)-pow(w6y,2.0));
  Phi(4,18) = -(dw6z+w6y*w6x)*sin(q6)-(pow(w6z,2.0)-pow(w6x,2.0))*cos(q6)-k65*(dw6x-w6z*w6y);
  Phi(4,19) = -(pow(w6y,2.0)-pow(w6z,2.0))*sin(q6)-(dw6z-w6y*w6x)*cos(q6)-k65*(dw6y+w6z*w6x);
  Phi(4,20) = w6z*(-w6y*sin(q6)+w6x*cos(q6))-k65*dw6z;
  Phi(4,21) = 0;
  Phi(4,22) = 0;
  Phi(4,23) = 0;
  Phi(4,24) = 0;
  Phi(4,25) = 0;
  Phi(4,26) = 0;
  Phi(4,27) = 0;
  Phi(4,28) = 0;
  Phi(4,29) = 0;
  Phi(4,30) = 0;
  Phi(4,31) = 0;
  Phi(4,32) = 0;
  Phi(4,33) = 0;
  Phi(4,34) = 0;
  Phi(4,35) = 0;
  Phi(4,36) = -cos(q5)*ddp4x-sin(q5)*ddp4y;
  Phi(4,37) = -w5x*w5z;
  Phi(4,38) = -dw5y;
  Phi(4,39) = -dw6x*sin(q6)-w6z*w6x*cos(q6)+k65*w6y*w6x;
  // Coppia: 6
  Phi(5,0)  = 0;
  Phi(5,1)  = 0;
  Phi(5,2)  = 0;
  Phi(5,3)  = 0;
  Phi(5,4)  = 0;
  Phi(5,5)  = 0;
  Phi(5,6)  = 0;
  Phi(5,7)  = 0;
  Phi(5,8)  = 0;
  Phi(5,9) = 0;
  Phi(5,10) = 0;
  Phi(5,11) = 0;
  Phi(5,12) = 0;
  Phi(5,13) = 0;
  Phi(5,14) = 0;
  Phi(5,15) = sin(q6)*ddp5x-cos(q6)*ddp5y;
  Phi(5,16) = cos(q6)*ddp5x+sin(q6)*ddp5y;
  Phi(5,17) = -pow(w6x,2.0)+pow(w6y,2.0);
  Phi(5,18) = -dw6x+w6z*w6y;
  Phi(5,19) = -dw6y-w6z*w6x;
  Phi(5,20) = -dw6z;
  Phi(5,21) = (k64-k54*k65)*ddq4-k65*ddq5-ddq6;
  Phi(5,22) = 0;
  Phi(5,23) = 0;
  Phi(5,24) = 0;
  Phi(5,25) = 0;
  Phi(5,26) = 0;
  Phi(5,27) = 0;
  Phi(5,28) = 0;
  Phi(5,29) = 0;
  Phi(5,30) = 0;
  Phi(5,31) = 0;
  Phi(5,32) = 0;
  Phi(5,33) = 0;
  Phi(5,34) = 0;
  Phi(5,35) = 0;
  Phi(5,36) = 0;
  Phi(5,37) = 0;
  Phi(5,38) = 0;
  Phi(5,39) = w6y*w6x;

  
  return Phi;
};
}
}
# endif
