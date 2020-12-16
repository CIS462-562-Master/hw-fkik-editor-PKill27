#include "aBehaviorController.h"

#include "aVector.h"
#include "aRotation.h"
#include <Windows.h>
#include <algorithm>



#define Truncate(a, b, c) (a = max<double>(min<double>(a,c),b))

double BehaviorController::gMaxSpeed = 200.0; //default = 1000.0, plugin = 6.0
double BehaviorController::gMaxAngularSpeed = 200.0; // default = 200.0 
double BehaviorController::gMaxForce = 2000.0;   // default = 2000.0
double BehaviorController::gMaxTorque = 2000.0;  // default = 2000.0
double BehaviorController::gKNeighborhood = 500.0; //default = 500.0, plugin = 100.0     
double BehaviorController::gOriKv = 32.0; //default = 32.0;    
double BehaviorController::gOriKp = 256.0; //default = 256.0;  
double BehaviorController::gVelKv = 10.0; //default = 10.0;    
double BehaviorController::gAgentRadius = 80.0; // default = 80.0, plugin = 2.0 


double BehaviorController::gMass = 1;
double BehaviorController::gInertia = 1;
double BehaviorController::KArrival = 1.0; 
double BehaviorController::KDeparture = 12000.0;
double BehaviorController::KNoise = 15.0;
double BehaviorController::KWander = 80.0;   
double BehaviorController::KAvoid = 600.0;  
double BehaviorController::TAvoid = 1000.0;   
double BehaviorController::KSeparation = 12000.0; 
double BehaviorController::KAlignment = 1.0;  
double BehaviorController::KCohesion = 1.0;  
double BehaviorController::leaderFollowDistance = 200.0;

const double M2_PI = M_PI * 2.0;

BehaviorController::BehaviorController() 
{
	m_state.resize(m_stateDim);
	m_stateDot.resize(m_stateDim);
	m_controlInput.resize(m_controlDim);

	vec3 m_Pos0 = vec3(0, 0, 0);
	vec3 m_Vel0 = vec3(0, 0, 0);
	vec3 m_lastVel0 = vec3(0, 0, 0);
	vec3 m_Euler = vec3(0, 0, 0);
	vec3 m_VelB = vec3(0, 0, 0);
	vec3 m_AVelB = vec3(0, 0, 0);
	
	m_Vdesired = vec3(0, 0, 0);
	m_lastThetad = 0.0;

	m_Active = true; 
	mpActiveBehavior = NULL;
	mLeader = false;
	m_leaderIndex = 0;

	reset();
}

AActor* BehaviorController::getActor()
{
	return m_pActor;
}

void BehaviorController::setActor(AActor* actor)

{
	m_pActor = actor;
	m_pSkeleton = m_pActor->getSkeleton();
}


void BehaviorController::createBehaviors(vector<AActor>& agentList, vector<Obstacle>& obstacleList)
{
	m_AgentList = &agentList;
	m_ObstacleList = &obstacleList;

	m_BehaviorList.clear();
	m_BehaviorList[SEEK] = new Seek(m_pBehaviorTarget);
	m_BehaviorList[FLEE] = new Flee(m_pBehaviorTarget);
	m_BehaviorList[ARRIVAL] = new Arrival(m_pBehaviorTarget);
	m_BehaviorList[DEPARTURE] = new Departure(m_pBehaviorTarget);
	m_BehaviorList[WANDER] = new Wander();
	m_BehaviorList[COHESION] = new Cohesion(m_AgentList);
	m_BehaviorList[ALIGNMENT] = new Alignment(m_pBehaviorTarget, m_AgentList);
	m_BehaviorList[SEPARATION] = new Separation(m_pBehaviorTarget, m_AgentList);
	m_BehaviorList[LEADER] = new Leader(m_pBehaviorTarget, m_AgentList);
	m_BehaviorList[FLOCKING] = new Flocking(m_pBehaviorTarget, m_AgentList);
	m_BehaviorList[AVOID] = new Avoid(m_pBehaviorTarget, m_ObstacleList);
}

BehaviorController::~BehaviorController()
{
	mpActiveBehavior = NULL;
}

void BehaviorController::reset()
{
	vec3 startPos;
	startPos[0] = ((double)rand()) / RAND_MAX;
	startPos[1] =  ((double)rand()) / RAND_MAX, 
	startPos[2] = ((double)rand()) / RAND_MAX;
	startPos = startPos - vec3(0.5, 0.5, 0.5);

	startPos[1] = 0.0; // set equal to zero for 2D case (assume y is up)
	mat3 Rmat;
	Rmat.Identity();  // sets initial orientation to be aligned with world coords

	m_Guide.setLocalTranslation(startPos * 500.0);
	m_Guide.setLocalRotation(Rmat);
	m_Guide.setGlobalTranslation(startPos * 500.0);
	m_Guide.setGlobalRotation(Rmat);
	
	for (int i = 0; i < m_stateDim; i++)
	{
		m_state[i] = 0.0;
		m_stateDot[i] = 0.0;
	}

	m_state[0] = m_Guide.getLocalTranslation();
	m_force = 0.0;
	m_torque = 0.0;
	m_thetad = 0.0;
	m_vd = 0.0;
}

///////////////////////////////////////////////////

inline void ClampAngle(double& angle)
{
	while (angle > M_PI)
	{
		angle -= M2_PI;
	}
	while (angle < -M_PI)
	{
		angle += M2_PI;
	}
}

void BehaviorController::sense(double deltaT)
{
	if (mpActiveBehavior)
	{
		// find the agents in the neighborhood of the current character.
	}
}

void BehaviorController::control(double deltaT)
// Given the active behavior this function calculates a desired velocity vector (Vdesired).  
// The desired velocity vector is then used to compute the desired speed (vd) and direction (thetad) commands
{
	if (mpActiveBehavior)
	{ 
		m_Vdesired = mpActiveBehavior->calcDesiredVel(this);
		m_Vdesired[1] = 0.0;

		//  force and torque inputs are computed from vd and thetad as follows:
		//              Velocity P controller : force = mass * Kv * (vd - v)
		//              Heading PD controller : torque = Inertia * (-Kv * thetaDot + Kp * (thetad - theta))
		//  where the values of the gains Kv and Kp are different for each controller
		// TODO: insert your code here

//#define OLDCONTROL

#ifdef OLDCONTROL
		//  Begin old implementation
		m_vd = m_Vdesired.Length();
		Truncate(m_vd, 0, gMaxSpeed);

		if (m_vd < 1.0)
		{
			m_vd = 0.0;
			m_thetad = m_lastThetad;
		}
		else
		{
			m_thetad = atan2(m_Vdesired[_X], m_Vdesired[_Z]);
			m_lastThetad = m_thetad;
		}
		ClampAngle(m_thetad);

		double theta = m_state[ORI][_Y];
		double deltaTheta = m_thetad - theta;
		ClampAngle(deltaTheta);

		m_force[0] = gMass * gVelKv * (0.0 - m_state[VEL][_X]);  // enforces no slip condition // 0.0;
		m_force[1] = 0.0;
		m_force[2] = gMass * gVelKv * (m_vd - m_state[VEL][_Z]);
		Truncate(m_force[2], -gMaxForce, gMaxForce);

		double Iyy = gInertia;
		m_torque[0] = 0.0;
		m_torque[1] = Iyy * (gOriKp * deltaTheta - gOriKv * m_state[AVEL][_Y]);
		m_torque[2] = 0.0;
		Truncate(m_torque[1], -gMaxTorque, gMaxTorque);

		if (m_vd < 2.0 &&  m_state[VEL][_Z] < 2.0)
		{
			m_force[2] = 0.0;
			m_torque[1] = 0.0;
		}
	}
	else
	{
		m_force[2] = 0.0;
		m_torque[1] = 0.0;
	}
	m_controlInput[0] = m_force;
	m_controlInput[1] = m_torque;

// end old implemenation

#else
// begin new implementation
		double speed_d;  // desired speed
		m_vd = speed_d = m_Vdesired.Length();
		Truncate(speed_d, 0, gMaxSpeed);

		if (speed_d < 1.0)
		{
			m_vd = speed_d = 0.0;
			m_thetad = m_lastThetad;
		}
		else
		{
			m_thetad = atan2(m_Vdesired[_X], m_Vdesired[_Z]);
			m_lastThetad = m_thetad;
		}
		ClampAngle(m_thetad);
		double thetad_degrees = m_thetad * (180.0 / M_PI);
	
	//////////////////////////////////////////////////////////////	
	// Force Controller

		vec3 w, w_x_V, Vb;

		// 3D Force Controller
		Vb = m_state[VEL];  // Vb is Body Axes Velocity
		w = m_state[AVEL];  // wb is Body Axes Angular Velocity
		w_x_V = w.Cross(Vb);

		m_force[0] = gMass * (w_x_V[0] + gVelKv * (0.0 - m_state[VEL][0]));  // no slip condition: Vb_x = 0
		m_force[1] = gMass * (w_x_V[1] + gVelKv * (0.0 - m_state[VEL][1]));  // want Vb_y = 0;
		m_force[2] = gMass * (w_x_V[2] + gVelKv * (speed_d - m_state[VEL][2]));   // speed_d is the desired speed
		Truncate(m_force[2], -gMaxForce, gMaxForce);

		// alternative force controller implementation using world velocity
		//mat3 Rmat;
		//Rmat.FromEulerAngles(mat3::YXZ, m_state[ORI]);
		//vec3 m_force = gMass * gVelKv * (m_Vdesired - Rmat* m_state[VEL]);

	//////////////////////////////////////////////////////////////	
	// Torque Controller
		mat3 I;
		// I = moment of inertia matrix (i.e. Inertia matrix) 
		// as default, assume I = Identity
		// If obtained from Unity Plugin I is in world coords. 
		// I = Rmat_trans * I * Rmat is conversion of I to body coords 

		I.Identity();
		I[1][1] = gInertia;
		double Iyy = I[1][1];
		double theta = m_state[ORI][_Y];
		double deltaTheta = m_thetad - m_state[ORI][_Y];
		ClampAngle(deltaTheta);

		// 2D Torque Controller implementation
		vec3 torque1;
		torque1[0] = 0.0;
		torque1[1] = (1/Iyy) * (gOriKp * deltaTheta - gOriKv * w[_Y]);
		torque1[2] = 0.0;
		Truncate(torque1[1], -gMaxTorque, gMaxTorque);

		/////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////
		//3D Controller Implementation

		// 1. Get current rotation matrix (Rmat)
		mat3 Rmat = m_Guide.getGlobalRotation(); 
		
		// 2. Construct desired rotation matrix (Rd)
		// Note: Z points forwards; Y Points UP; and X points left
		vec3 Xaxis, Yaxis, Zaxis;

		// compute Z column of Rd from direction of desired velocity vector
		if (speed_d > 0.0001) // handles case where Vdesired is close to zero
		{
			Zaxis = m_Vdesired;
			Zaxis.Normalize();
		} 
		else Zaxis = vec3 (Rmat[0][2], Rmat[1][2], Rmat[2][2]);

		// compute Y column of Rd given Up = vec3(0, 1, 0)
		vec3 up(0.0, 1.0, 0.0);
		Yaxis = up;

		// compute X column of Rd given Zaxis and Yaxis vectors
		Xaxis = Yaxis.Cross(Zaxis);

		mat3 Rd(Xaxis, Yaxis, Zaxis);
		Rd = Rd.Transpose();  //mat3 constructor above sets rows of Rd matrix.  Transpose to put in column format

		// compute the Rd euler angles
		vec3 theta_d;
		Rd.ToEulerAngles(mat3::YXZ, theta_d);
		for (int i = 0; i< 3; i++)
			ClampAngle(theta_d[i]);
		
		// 3. Compute the change in rotation (delR) that will achieve Rd  (Rd = Rmat * delR)
		vec3 delTheta, Theta;
		Rmat.ToEulerAngles(mat3::YXZ, Theta); 

		mat3 delR = Rmat.Transpose() * Rd;
		delR.Reorthogonalize();
		delR.ToEulerAngles(mat3::YXZ, delTheta);
		for (int i=0; i< 3; i++)
			ClampAngle(delTheta[i]);

		// 4. get corresponding Axis/Angle representation associated with delR
		vec3 axis;   
		double angle;
		delR.ToAxisAngle(axis, angle);   

		if (abs(angle) < 0.0001)
		{
			angle = 0.0;
			axis = vec3(0.0, 1.0, 0.0);
		}

		// 5. compute torque = I * (Kp*axis*angle - Kv *w) + w x Iw
		vec3 w_x_Iw, Iw;
		Iw = I * w;
		w_x_Iw = w.Cross(Iw);

		vec3 torque2 = I * (gOriKp * axis * angle - gOriKv * w) + w_x_Iw;
		Truncate(torque2[1], -gMaxTorque, gMaxTorque);

		m_torque = torque2;
		assert(!isnan(m_torque[0]));  assert(!isnan(m_torque[1]));  assert(!isnan(m_torque[2]));

		// set forces and torques to zero for if essentially not moving
		if (speed_d < 2.0 &&  Vb.Length() < 2.0)
		{
			m_force[2] = 0.0;
			m_torque[1] = 0.0;
		}
	}
	else // set forces and torques to zero if behavior not valid
	{
		m_force[2] = 0.0;
		m_torque[1] = 0.0;
	}
	m_controlInput[0] = m_force;
	m_controlInput[1] = m_torque;
#endif
}

void BehaviorController::act(double deltaT)
{
	computeDynamics(m_state, m_controlInput, m_stateDot, deltaT);
	
	int EULER = 0;
	int RK2 = 1;
	updateState(deltaT, EULER);
}


void BehaviorController::computeDynamics(vector<vec3>& state, vector<vec3>& controlInput, vector<vec3>& stateDot, double deltaT)
// Compute stateDot vector given the control input and state vectors
//  This function sets derive vector to appropriate values after being called
{
	vec3& force = controlInput[0];
	vec3& torque = controlInput[1];

//#define OLDDYNAMICS
#ifdef  OLDDYNAMICS
// old implementation
	stateDot[POS][_X] = state[VEL][_Z] * sin(m_state[ORI][_Y]);
	stateDot[POS][_Y] = 0.0;
	stateDot[POS][_Z] = state[VEL][_Z] * cos(m_state[ORI][_Y]);
	stateDot[ORI] = state[AVEL];
	stateDot[VEL] = force / gMass;
	stateDot[AVEL] = torque / gInertia;
#else
// begin new Dynamics Implementation
	vec3 w_x_V, Vb, w, w_x_Iw, Iw, theta;

	mat3 Rmat = m_Guide.getGlobalRotation();

	// I = moment of inertia matrix (i.e. Inertia matrix) 
	// as default, assume I = Identity
	// If obtained from Unity Plugin I is in world coords. 
	// I = Rmat_trans * I * Rmat is conversion of I to body coords 
	mat3 I, I_inv;  I.Identity(); I_inv = I.Inverse(); 

	Vb = state[VEL];
	w = state[AVEL];  // wb is Body Axes Angular Velocity
	w_x_V = w.Cross(Vb);
	Iw = I * w;
	w_x_Iw = w.Cross(Iw);

	double theta_x, theta_y, theta_z;
	theta = m_state[ORI];

	mat3 L, Linv;
	L[0][0] = cos(theta[2]);  L[0][1] = cos(theta[0])*sin(theta[2]);  L[0][2] = 0.0;
	L[1][0] = -sin(theta[2]); L[1][1] = cos(theta[0])*cos(theta[2]);  L[1][2] = 0.0;
	L[2][0] = 0.0;            L[2][1] = -sin(theta[0]);               L[2][2] = 1.0;

	Linv = L.Inverse();

	stateDot[POS] = Rmat * Vb;                   // This transform body velocity to world velocity
	stateDot[ORI] = Linv * w;                    // angular velocity w is in Body Axes
	stateDot[VEL] = (1/ gMass) * force - w_x_V;  //
	stateDot[AVEL] = I_inv*(torque - w_x_Iw);    //

//end new implementation

#endif
}

void BehaviorController::updateState(float deltaT, int integratorType)
//  Update the state vector given the m_stateDot vector
//  Perform validation check to make sure all values are within MAX values
{
	switch (integratorType)
	{
		case 0:
			for (int i = 0; i < m_stateDim; i++)
				m_state[i] += m_stateDot[i] * deltaT;
			break;

		case 1:
		{
			vector<vec3> predictedState;
			predictedState.resize(m_stateDim);

			vector<vec3> predictedStateDot;
			predictedStateDot.resize(m_stateDim);

			predictedState = m_state;
			for (int i = 0; i < m_stateDim; i++)
				predictedState[i] += m_stateDot[i] * deltaT;

			computeDynamics(predictedState, m_controlInput, predictedStateDot, deltaT);

			for (int i = 0; i < m_stateDim; i++)
			{
				m_stateDot[i] = 0.5*(m_stateDot[i] + predictedStateDot[i]);
				m_state[i] += m_stateDot[i] * deltaT;
			}

			break;
		}
	}

//#define OLDUPDATE

#ifdef OLDUPDATE
	// Begin old method to update Guide Orientation (rotation matrix = Rmat)

	// compute direction from nonzero velocity vector.  need to handle case where dir = 0.0;
	vec3 dir;
	double speed = m_Vel0.Length();
	if (speed < 0.00001)
	dir = vec3(0.0, 0.0, 1.0);
	else if (m_Vel0.Length() < 1.0)
	{
	dir = m_lastVel0;
	dir.Normalize();;
	m_state[ORI][_Y] = atan2(dir[_X], dir[_Z]);
	}
	else
	{
	dir = m_Vel0;
	m_lastVel0 = m_Vel0;
	}

	// Compute current rotation matrix R from knowledge of world velocity  dir = m_Vel0
	dir.Normalize();  // forward (z) vector with respect to world coordinates
	vec3 up(0.0, 1.0, 0.0);
	vec3 left = up.Cross(dir);
	left.Normalize();

	mat3 R(left, up, dir);  //this sets rows of rotation R.  not standard column format
	R = R.Transpose();   // columns of R now represent x, y and z coordinate axes
	m_Guide.setLocalRotation(Rmat);
	m_Guide.setGlobalRotation(Rmat);
	m_Guide.setLocalTranslation(m_Guide.getLocalTranslation() + m_Vel0 * deltaT);
	m_Guide.setGlobalTranslation(m_Guide.getGlobalTranslation() + m_Vel0 * deltaT);
	*/
	// end old implementation

#else

	m_Pos0 = m_state[POS];
	m_Euler = m_state[ORI];
	m_VelB = m_state[VEL];
	m_AVelB = m_state[AVEL];
	m_Vel0 = m_stateDot[POS];
	for (int i = 0; i < 3; i++)
	{
		ClampAngle(m_state[ORI][i]);
		Truncate(m_state[VEL][i], -gMaxSpeed, gMaxSpeed);
		Truncate(m_state[AVEL][i], -gMaxAngularSpeed, gMaxAngularSpeed);
	}
	mat3 Rmat;
	Rmat.FromEulerAngles(mat3::YXZ, m_Euler);
	m_Guide.setLocalRotation(Rmat);
	m_Guide.setGlobalRotation(Rmat);
	m_Guide.setLocalTranslation(m_Guide.getLocalTranslation() + m_Vel0 * deltaT);
	m_Guide.setGlobalTranslation(m_Guide.getGlobalTranslation() + m_Vel0 * deltaT);

#endif
	
}


void BehaviorController::setTarget(AJoint& target)
{
	m_pBehaviorTarget = &target;
	for (unsigned int i = 0; i < m_BehaviorList.size(); i++)
	{
		BehaviorType index = (BehaviorType) i;
		m_BehaviorList[index]->setTarget(m_pBehaviorTarget);
	}
}

void BehaviorController::setActiveBehavior(Behavior* pBehavior)
{
	mpActiveBehavior = pBehavior;
}

void BehaviorController::setActiveBehaviorType(BehaviorType type)
{
	m_BehaviorType = type;
	Behavior* pActiveBehavior = m_BehaviorList[type];
	setActiveBehavior(pActiveBehavior);

}

vec3 BehaviorController::getForce(RefFrameType type)
{
	vec3 force;
	if (type == WORLD)
		 force = m_Guide.getGlobalRotation()*m_force;
	else force = m_force;

	return force;
}

vec3 BehaviorController::getTorque(RefFrameType type)
{
	vec3 torque;
	if (type == WORLD)
		torque = m_Guide.getGlobalRotation()*m_torque;
	else torque = m_torque;

	return torque;
}

void BehaviorController::display()
{ // helps with debugging behaviors.  red line is actual velocity vector, green line is desired velocity vector
#ifdef DEBUG_BEHAVIORDISPLAY	
	vec3 pos = getPosition();
	double scale = 1.0;
	vec3 vel = scale* getVelocity();
	double velMag = vel.Length();
	vec3 dvel = scale* getDesiredVelocity();
	vec3 angle = getOrientation() * (180.0 / 3.14159);

	glBegin(GL_LINES);
	glColor3f(1, 0, 0);
	glVertex3f(pos[0], pos[1], pos[2]);
	glVertex3f(pos[0] + vel[0], pos[1] + vel[1], pos[2] + vel[2]);
	glColor3f(0, 1, 0);
	glVertex3f(pos[0], pos[1], pos[2]);
	glVertex3f(pos[0] + dvel[0], pos[1] + dvel[1], pos[2] + dvel[2]);
	glEnd();

	if (this->isLeader())
		glColor3f(0, 0, 1);
	else glColor3f(0.5, 0, 0);

	glPushMatrix();
	glTranslatef(pos[0], pos[1], pos[2]);
	//glRotatef(90 - angle[1], 0, 1, 0);
	glRotatef(angle[1], 0, 1, 0);
	glutSolidCone(40, 80, 10, 10);
	glutSolidSphere(35, 10, 10);
	glPopMatrix();

	BehaviorType active = getActiveBehaviorType();
	Behavior* pBehavior = m_BehaviorList[active];
	pBehavior->display(this);
#endif
}


