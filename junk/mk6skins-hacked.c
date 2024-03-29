#include "posemath.h"         /* save as mk6skins.c */
#include "rtapi_math.h"
#include "kinematics.h"       /* decls for kinematicsForward, etc. */
#include "sincos.h"
#include "rtapi.h"		/* RTAPI realtime OS API */
#include "rtapi_app.h"		/* RTAPI realtime module decls */
#include "hal.h"

double abanglefromtriangleabc(double a, double b, double c)
{
    return acos((a*a + b*b - c*c)/(2*a*b)); 
}

double clengthfromtriangleabangle(double a, double b, double ang)
{
    return sqrt(a*a + b*b - 2*a*b*cos(ang)); 
}

struct arcdata_data {
    hal_float_t a, b, c, d, e, f, ah, tblx, armx;
    hal_float ac, ab; 
};

void setupconstants(struct arcdata_data *hd) 
{
    hd->a = 743.5035783400; 
    hd->b = 754.8024952389; 
    hd->c = 745.4991150306; 
    hd->d = 697.7062405262; 
    hd->e = 695.5193939239; 
    hd->f = 583.6487506837; 
    hd->ah = 0.0764480526; 
    hd->tblx = 525.4558306022; 
    hd->armx = 396.1544040232; 
}

void setupprecalcs(struct arcdata_data *hd) 
{
    hd->ab = abanglefromtriangleabc(a, b, c); 
    hd->ac = abanglefromtriangleabc(a, c, b); 
}

void jointstoxy(struct arcdata_data* hd, const double* joint, EmcPose* world)
{
    double tbl, arm;
    double aARM, aTBL, avec, pvec;

    tbl = joint[0] + hd->tblx;
    arm = joint[1] + hd->armx;
    
    aTBL = abanglefromtriangleabc(hd->b, hd->d, tbl); 
    aARM = abanglefromtriangleabc(hd->c, hd->e, arm); 
    avec = hd->ab - aTBL; 
    pvec = avec + hd->ac - aARM + hd->ah; 

    world->tran.x =     0 - hd->a*sin(avec) + hd->f*sin(pvec); 
    world->tran.y = hd->d - hd->a*cos(avec) + hd->f*cos(pvec); 
    world->tran.z = joint[2];
}

void xytojoints(struct arcdata_data* hd, const EmcPose* world, double* joint)
{
    double tbl, arm;    
    double gx, gy, g; 
    double aTBL, aARM, af, avec, pvec;

    gx = world->tran.x;
    gy = world->tran.y - hd->d;
    g = sqrt(gx*gx + gy*gy);  // line from corner ab to p

    // consider a vertical line through corner ab, splitting angle ag and parallel to vertical line through corner af. use parallel axiom
    avec = abanglefromtriangleabc(hd->a, g, hd->f) - atan2(gx, -gy); 
    af = abanglefromtriangleabc(hd->a, hd->f, g); 
    pvec = avec + af; 
    aTBL = hd->ab - avec; 
    aARM = hd->ac - af + hd->ah; 

    tbl = clengthfromtriangleabangle(hd->d, hd->b, aTBL); 
    arm = clengthfromtriangleabangle(hd->c, hd->e, aARM); 

    joint[0] = tbl - hd->tblx;
    joint[1] = arm - hd->armx;
    joint[2] = world->tran.z;
}

/* Python translation

from math import acos, cos, sin, atan2, sqrt

def abanglefromtriangleabc(a, b, c):
    return acos((a*a + b*b - c*c)/(2*a*b)) 
def clengthfromtriangleabangle(a, b, ang):
    return sqrt(a*a + b*b - 2*a*b*cos(ang)) 

class HalData:
    def __init__(self):
        self.a = 743.5035783400; 
        self.b = 754.8024952389; 
        self.c = 745.4991150306; 
        self.d = 697.7062405262; 
        self.e = 695.5193939239; 
        self.f = 583.6487506837; 
        self.ah = 0.0764480526; 
        self.tblx = 525.4558306022; 
        self.armx = 396.1544040232; 

        self.ab = abanglefromtriangleabc(self.a, self.b, self.c) 
        self.ac = abanglefromtriangleabc(self.a, self.c, self.b) 

    def jointstoxy(self, joint):
        tbl = joint[0] + self.tblx
        arm = joint[1] + self.armx

        aTBL = abanglefromtriangleabc(self.b, self.d, tbl)
        aARM = abanglefromtriangleabc(self.c, self.e, arm) 
        
        avec = self.ab - aTBL; 
        pvec = avec + self.ac - aARM + self.ah 

        return (0      - self.a*sin(avec) + self.f*sin(pvec),  
                self.d - self.a*cos(avec) + self.f*cos(pvec))

    def xytojoints(self, tran):
        gx = tran[0]
        gy = tran[1] - self.d
        g = sqrt(gx*gx + gy*gy)
        
        avec = abanglefromtriangleabc(self.a, g, self.f) - atan2(gx, -gy) 
        af = abanglefromtriangleabc(self.a, self.f, g) 
        pvec = avec + af 
        aTBL = self.ab - avec 
        aARM = self.ac - af + self.ah
        
        tbl = clengthfromtriangleabangle(self.b, self.d, aTBL) 
        arm = clengthfromtriangleabangle(self.c, self.e, aARM) 

        return (tbl - self.tblx, arm - self.armx)
*/

struct arcdata_data *haldata; // global object with all the parameters in it

int kinematicsForward(const double* joint, EmcPose* world, const KINEMATICS_FORWARD_FLAGS* fflags, KINEMATICS_INVERSE_FLAGS* iflags)
{
    jointstoxy(haldata, joint, world); 
    return (0);
}

int kinematicsInverse(const EmcPose* world, double* joint, const KINEMATICS_INVERSE_FLAGS* iflags, KINEMATICS_FORWARD_FLAGS* fflags)
{
    xytojoints(haldata, world, joint); 
    *fflags = 0;
    return (0);
}

int kinematicsHome(EmcPose* world, double*joint, KINEMATICS_FORWARD_FLAGS* fflags, KINEMATICS_INVERSE_FLAGS* iflags)
{
    *fflags = 0;
    *iflags = 0;
    return kinematicsForward(joint, world, fflags, iflags);
}

KINEMATICS_TYPE kinematicsType()
{
    return KINEMATICS_BOTH;
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
    setupconstants(haldata); 
    setupprecalcs(haldata); // should be called in case a value is changed from the hal_param system below

    if (hal_param_float_new("kins-a", HAL_RW, &(haldata->a), comp_id) < 0){
        rtapi_print("failed to make pin");
        return -EINVAL;
    }
    if (hal_param_float_new("kins-b", HAL_RW, &(haldata->b), comp_id) < 0){
        rtapi_print("failed to make pin");
        return -EINVAL;
    }
    if (hal_param_float_new("kins-c", HAL_RW, &(haldata->c), comp_id) < 0){
        rtapi_print("failed to make pin");
        return -EINVAL;
    }
    if (hal_param_float_new("kins-d", HAL_RW, &(haldata->d), comp_id) < 0){
        rtapi_print("failed to make pin");
        return -EINVAL;
    }
    if (hal_param_float_new("kins-e", HAL_RW, &(haldata->e), comp_id) < 0){
        rtapi_print("failed to make pin");
        return -EINVAL;
    }
    if (hal_param_float_new("kins-f", HAL_RW, &(haldata->f), comp_id) < 0){
        rtapi_print("failed to make pin");
        return -EINVAL;
    }
    if (hal_param_float_new("kins-ab", HAL_RW, &(haldata->ab), comp_id) < 0){
        rtapi_print("failed to make pin");
        return -EINVAL;
    }
    if (hal_param_float_new("kins-ac", HAL_RW, &(haldata->ac), comp_id) < 0){
        rtapi_print("failed to make pin");
        return -EINVAL;
    }
    if (hal_param_float_new("kins-ah", HAL_RW, &(haldata->ah), comp_id) < 0){
        rtapi_print("failed to make pin");
        return -EINVAL;
    }
    if (hal_param_float_new("kins-tblx", HAL_RW, &(haldata->tblx), comp_id) < 0){
        rtapi_print("failed to make pin");
        return -EINVAL;
    }
    if (hal_param_float_new("kins-armx", HAL_RW, &(haldata->armx), comp_id) < 0){
        rtapi_print("failed to make pin");
        return -EINVAL;
    }
    hal_ready(comp_id);
    return 0;
}

void rtapi_app_exit(void) { hal_exit(comp_id); }



