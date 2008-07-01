#include <iostream>

using namespace std;

/*#include <DTRA/fc.h>*/

/*  FORTRAN pseudo-code

    parameter MXND = 100000
    double precision crd, vel, acc, frc, dt
    dimensions crd(3*MXND), vel(3*MXND), acc(3*MXND), frc(3,*MXND)
c   -- I believed this should work too
    dimensions crd(3,MXND), vel(3,MXND), acc(3,MXND), frc(3, MXND)
c   Initialization:
    integer opt
    opt = 1
    mxn = MXND
    dt = 0.;
c   -- mars returns  nnd and crd
    call mars3D_RPICFD(opt, dt, mxn, nnd, crd, vel, acc, frc)
c   -- loop
    dt = ...
    while (true) 
       if (.not.converged) 
          opt = 16
c         -- store velocities in vel and accelerations in acc
          call mars3D_RPICFD(opt, dt, mxn, nnd, crd, vel, acc, frc)
c         -- retrieve forces from frc
       else
          opt = 32
c         -- store velocities in vel and accelerations in acc
          call mars3D_RPICFD(opt, dt, mxn, nnd, crd, vel, acc, frc)
c         -- retrieve forces from frc
       endif
    end
c   -- termination (not necessary)
    stp = 12
    call mars3D_RPICFD(opt, dt, mxn, nnd, crd, vel, acc, frc)
*/

extern "C" void mars3d_RPICFD (
        int *i4opt,       // options 
                          // =1 read input file,
                          //    initialize,
                          //    return node locations
                          // =2 read restart file,
                          //    initialize,
                          //    return node locations
                          // =4 write restart file,
                          // =8 stop execution
                          // =12 write restart file and stop execution
                          // =16 reset state variable at beginning of step
                          //     compute nodal forces
                          //     return forces
                          // =32 save state variables
                          //     advance one step
                          //     compute nodal forces
                          //     return forces
        double *dt,       // time interval
        int *m4wn,        // maximum number of nodal points (from CFD)
        int *n4wn,        // current number of nodal points
        double *c8wn,     // node coordinates (to CFD first time only)
        double *v8wn,     // node velocities (from CFD every step)
        double *a8wn,     // node acceleractions (from CFD every step)
        double *f8wn      // node forces (to CFD every step)
        ) {
    int opt = *i4opt;
    if (opt&1) {
        cout << "Read and initialize problem " << endl;
        *n4wn = 10;
        if (*n4wn > *m4wn) {
            cout << "Increase array size" << endl;
            return;
        }
        for (int i=0; i<30; i++) 
            c8wn[i] = (double)(i+1.);
        // c1x=1., c1y=2. c1z=3. c2x=4. c2y=5. c2z=6. c3x=7. etc
    } else if (opt&16 || opt&32) {
        cout << "Receive velocities and accelerations" << endl;
        cout << "Velocities" << endl;
        for (int i=0, j=0; i<10; i++) 
            cout << "Node " << i << "(cx, cy, cz): "
                 << v8wn[j++] << ", " << v8wn[j++] << ", " 
                 << v8wn[j++] << endl;
        cout << "Accelerations" << endl;
        for (int i=0, j=0; i<10; i++) 
            cout << "Node " << i << "(vx, vy, vz): "
                 << a8wn[j++] << ", " << a8wn[j++] << ", " 
                 << a8wn[j++] << endl;
        for (int i=0; i<30; i++) 
            f8wn[i] = (double)(i+1.);
        // f1x=1., f1y=2. f1z=3. f2x=4. f2y=5. f2z=6. f3x=7. etc
    } else 
        cout << "Invalid option" << endl;
    return;
}
