#include "posemath.h"         /* save as mk6skins.c */
#include "rtapi_math.h"
#include "kinematics.h"       /* decls for kinematicsForward, etc. */
#include "sincos.h"
#include "rtapi.h"		/* RTAPI realtime OS API */
#include "rtapi_app.h"		/* RTAPI realtime module decls */
#include "hal.h"

struct arcdata_data {
    hal_float_t a, b, c, d, e, f, ac, ab;   /* double constants */
} *haldata;

/* constants*/

#define A  (haldata->a)
#define B  (haldata->b)
#define C  (haldata->c)
#define D  (haldata->d)
#define E  (haldata->e)
#define F  (haldata->f)
#define AC (haldata->ac)
#define AB (haldata->ab)

/* joint[0],[1] and [2] linear units */

int kinematicsForward(const double * joint,
                      EmcPose * world,
                      const KINEMATICS_FORWARD_FLAGS * fflags,
                      KINEMATICS_INVERSE_FLAGS * iflags)
{
    double tbl, arm;                       /* variables*/
    double aARM, aX, aF, aG, aTBL;
    double x, y, g;
/*in equations use A,B,C,D,E,F,AB,AC,x,y,g,aX,aF,aG,aTBL,aARM,tbl,arm */

    tbl = joint[0];                                   /*<<<<<<< tbl*/
    arm = joint[1];                                   /*<<<<<<< arm*/
    
    aTBL = acos((((D*D)+(B*B))-(tbl*tbl))/(2*D*B));
    aARM = acos((((C*C)+(E*E))-(arm*arm))/(2*C*E)); 
    aG   = AB-aARM+0.0799634827644152;
    aF   = atan((F*sin(aG))/(A-(F*cos(aG))));
    aX   = aTBL-(AC-aF)-0.4837;
    g    = sqrt (((F*F)+(A*A))-((2*F*A)*( cos(aG))));
    x    = (g*sin(aX))+110;
    y    = 572-(g*cos(aX));
	


    world->tran.x = x;                                /*>>>>>x*/
    world->tran.y = y;                                /*>>>>>y*/
    world->tran.z = joint[2];                         /*z>>>>>*/
  
    return (0);
}


int kinematicsInverse(const EmcPose * world,
                      double * joint,
                      const KINEMATICS_INVERSE_FLAGS * iflags,
                      KINEMATICS_FORWARD_FLAGS * fflags)
{
    double tbl, arm;    /* variables*/
    double x, y;
    double g, aX, aF, aG, aTBL, aARM;
/*in equations use A,B,C,D,E,F,AB,AC,x,y,g,aX,aF,aG,aTBL,aARM,tbl,arm */

    x = world->tran.x;                                         /*'<<<<<<x*/
    y = world->tran.y;                                         /*'<<<<<<y*/
    joint[2] = world->tran.z;                                  /*'z>>>>>>*/
 
    g = sqrt(((572-y)*(572-y))+((x-110)* (x-110)));
    aX = asin((x-110)/g); /*check for asin command???????*/
    aF = acos((((A*A)+(g*g))-(F *F ))/(2*A*g));
    aTBL = AC-aF+aX+0.4837;
    tbl = sqrt(((D*D)+(B*B))-((2*D*B)*(cos(aTBL))));
    aG = acos((((F*F)+(A*A))-( g*g ))/(2*F*A));
    aARM = AB-aG+0.0799634827644152;
    arm = sqrt(((C*C)+(E*E))-((2*C*E)*(cos(aARM))));



    joint[0] = tbl;                                  /*'tbl>>>>>>*/
    joint[1] = arm;                                  /*'arm>>>>>>*/
    
    *fflags = 0;

    return (0);
}


/* Kinematics Both */

/* constants*/
#define DEFAULT_A 744.942619
#define DEFAULT_B 746.50577
#define DEFAULT_C 743.340584
#define DEFAULT_D 695.063413
#define DEFAULT_E 692.996544
#define DEFAULT_F 582.632544
#define DEFAULT_AB 1.05086603051362;  /* 60degrees */
#define DEFAULT_AC 1.04350675226089;

    

/* implemented for these kinematics as giving joints preference */
int kinematicsHome(EmcPose * world,
		   double *joint,
		   KINEMATICS_FORWARD_FLAGS * fflags,
		   KINEMATICS_INVERSE_FLAGS * iflags)
{
    *fflags = 0;
    *iflags = 0;

    return kinematicsForward(joint, world, fflags, iflags);
}

KINEMATICS_TYPE kinematicsType()
{    return KINEMATICS_BOTH;
}


EXPORT_SYMBOL(kinematicsType);
EXPORT_SYMBOL(kinematicsForward);
EXPORT_SYMBOL(kinematicsInverse);

int comp_id;
int rtapi_app_main(void) {
    comp_id = hal_init("mk6skins");
    if(comp_id < 0) {
        return comp_id;
    }
    haldata = hal_malloc(sizeof(struct arcdata_data));
    if (hal_param_float_new("kins-a", HAL_RW, &(haldata->a), comp_id) < 0){
        rtapi_print("failed to make pin");
        return -EINVAL;
    }
    A = DEFAULT_A;
    if (hal_param_float_new("kins-b", HAL_RW, &(haldata->b), comp_id) < 0){
        rtapi_print("failed to make pin");
        return -EINVAL;
    }
    B = DEFAULT_B;
    if (hal_param_float_new("kins-c", HAL_RW, &(haldata->c), comp_id) < 0){
        rtapi_print("failed to make pin");
        return -EINVAL;
    }
    C = DEFAULT_C;
    if (hal_param_float_new("kins-d", HAL_RW, &(haldata->d), comp_id) < 0){
        rtapi_print("failed to make pin");
        return -EINVAL;
    }
    D = DEFAULT_D;
    if (hal_param_float_new("kins-e", HAL_RW, &(haldata->e), comp_id) < 0){
        rtapi_print("failed to make pin");
        return -EINVAL;
    }
    E = DEFAULT_E;
    if (hal_param_float_new("kins-f", HAL_RW, &(haldata->f), comp_id) < 0){
        rtapi_print("failed to make pin");
        return -EINVAL;
    }
    F = DEFAULT_F;
    if (hal_param_float_new("kins-ab", HAL_RW, &(haldata->ab), comp_id) < 0){
        rtapi_print("failed to make pin");
        return -EINVAL;
    }
    AB = DEFAULT_AB;
    if (hal_param_float_new("kins-ac", HAL_RW, &(haldata->ac), comp_id) < 0){
        rtapi_print("failed to make pin");
        return -EINVAL;
    }
    AC = DEFAULT_AC;
    hal_ready(comp_id);
    return 0;
}

void rtapi_app_exit(void) { hal_exit(comp_id); }



