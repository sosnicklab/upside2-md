#include "deriv_engine.h"
#include <string>
#include "timing.h"
#include "affine.h"
#include <cmath>
#include "h5_support.h"
#include <vector>
#include "spline.h"
#include <iostream>
#include <map>
#include <set>

#define N_AMINO_ACID_TYPE 20
#define ForceNorm 0.70710678
#define FORCE_UNIT 2.428602e-5

using namespace h5;
using namespace std;

typedef vector<int>   int1d;
typedef vector<int1d> int2d;
typedef vector<int2d> int3d;
typedef map<int, int1d> map_int1d;

void RotateXY(vector<float>& x, vector<float>& y, float rx, float ry, float angle, int size) {
    for (int na=0; na< size; ++na ) {
        float dx = x[na] - rx;
        float dy = y[na] - ry;
        if (dx==0 and dy==0) continue;
        x[na]    = rx + cos(angle) * dx - sin(angle) * dy;
        y[na]    = ry + sin(angle) * dx + cos(angle) * dy;
    }
}

void FindMinMax(vector<float>& X, float& minx, float& maxx, bool initialization, int size ) {
    int ii = 0;
    for (float& x: X) {
        if (ii >= size) break;
        if(not initialization) {
            minx = x;
            maxx = x;
            initialization = true;
        }
        else {
            if(x > maxx)
                maxx = x;
            else if (x< minx)
                minx = x;
        }
        ii++;
    }
}

struct Surface : public CoordNode
{
    struct auxiliary_point { float xz[2]; };

    int method_type; // 0: fast algorithm; 1: diffusion algorithm

    int n_residue;  // number of residues
    int n_used_res;
    int n_rotation; // The number of rotations, from 0 to 90 degrees
    int n_aux;     
    float r_n_rotation; // The number of rotations, from 0 to 90 degrees
    float dangle;   // the rotation angle of every rotation

    float half_thickness;
    float zbuffer;
    float atom_radius;
    float tm_min;
    float tm_max;

    int NGRIDZ;
    float GRID_X;
    float GRID_Y;
    float GRID_Z;
    float r_GRID_Z;
    int SNGRID_X;
    int SNGRID_Y;
    int SNGRID_Z;

    CoordNode& pos;  // N, CA, CB C atoms

    int1d index_cb;
    int1d type_count;
    int1d used_list;
    int1d tm_list;

    vector<float> cbx;
    vector<float> cby;
    vector<float> cbz;
    vector<float> bbx;
    vector<float> bby;
    vector<float> bbz;
    vector<float> fx;
    vector<float> fy;
    vector<float> sfx;
    vector<float> sfy;

    vector<auxiliary_point> auxiliary_points;

    //read_attribute<int>(grp, ".", "n_dim1")+read_attribute<int>(grp, ".", "n_dim2")),
    Surface(hid_t grp, CoordNode& pos_):
        CoordNode( get_dset_size(1,grp,"index_cb")[0], 3 ),

        method_type( read_attribute<int>(grp, ".", "method_type") ),
        n_residue( get_dset_size(1,grp,"index_cb")[0] ),
        n_used_res(n_residue),
        n_rotation( read_attribute<int>(grp, ".", "n_rotation") ),
        n_aux(get_dset_size(2,grp,"auxiliary_point")[0]),
        r_n_rotation( 1./n_rotation ),
        dangle(0.5*M_PI_F*r_n_rotation),
        half_thickness(read_attribute<float>(grp, ".", "half_thickness")),
        zbuffer(read_attribute<float>(grp, ".", "zbuffer")),
        atom_radius( read_attribute<float>(grp, ".", "atom_radius") ),
        tm_min(-half_thickness-zbuffer),
        tm_max(half_thickness+zbuffer),

        NGRIDZ( read_attribute<int>(grp, ".", "n_gridz") ),
        GRID_X(read_attribute<float>(grp, ".", "len_gridx")),
        GRID_Y(read_attribute<float>(grp, ".", "len_gridy")),
        GRID_Z((tm_max-tm_min+2*atom_radius+0.1)/NGRIDZ),
        r_GRID_Z(1.f/GRID_Z),

        SNGRID_X( read_attribute<int>(grp, ".", "n_smallgridx") ),
        SNGRID_Y( read_attribute<int>(grp, ".", "n_smallgridy") ),
        SNGRID_Z( read_attribute<int>(grp, ".", "n_smallgridz") ),

        pos(pos_),
        index_cb(n_residue),
        type_count(n_residue, 0),
        used_list(n_residue),
        tm_list(n_residue),

        cbx(n_residue),
        cby(n_residue),
        cbz(n_residue),
        bbx(n_residue*3),
        bby(n_residue*3),
        bbz(n_residue*3),
        fx(n_residue,  0.),
        fy(n_residue,  0.),
        sfx(n_residue, 0.),
        sfy(n_residue, 0.),

        auxiliary_points(n_aux)

    {
        check_elem_width_lower_bound(pos, 3);
        traverse_dset<1,   int>(grp, "index_cb",        [&](size_t rt, index_t x)      { index_cb[rt] = x;});
        traverse_dset<2, float>(grp, "auxiliary_point", [&](size_t rt, int d, float x) { auxiliary_points[rt].xz[d] = x;});

        // 0, no exclusion; 
        // 1, exclusion by exclusion_list; 
        int exclude_way = read_attribute<int>(grp, ".", "exclusion"); 

        if (exclude_way == 0) {
            for(int nr: range(n_residue)) used_list[nr] = nr;
        }
        else if  (exclude_way == 1) {
            n_used_res = get_dset_size(1,grp,"included_list")[0];
            used_list.resize(n_used_res);
            tm_list.resize(n_used_res);
            traverse_dset<1, int>(grp, "included_list", [&](size_t ne, index_t rt){ used_list[ne] = rt;});
        }

        if (n_used_res != n_residue) {
            cbx.resize(n_used_res);
            cby.resize(n_used_res);
            cbz.resize(n_used_res);
            bbx.resize(n_used_res*3);
            bby.resize(n_used_res*3);
            bbz.resize(n_used_res*3);
            fx.resize(n_used_res);
            fy.resize(n_used_res);
        }
    }

    virtual void compute_value(ComputeMode mode) {

        Timer timer(string("surface"));
        VecArray all_pos = pos.output;

        // 0. find the TM residues
        int ii = 0;
        for(int nr: used_list) {
            float z1 = all_pos(2, nr*3+0);
            if (z1 < tm_max and z1 > tm_min) {
                float z2 = all_pos(2, nr*3+1);
                if (z2 < tm_max and z2 > tm_min) {
                    float z3 = all_pos(2, nr*3+2);
                    if (z3 < tm_max and z3 > tm_min) {
                        float z4 = all_pos(2, index_cb[nr]);
                        if (z4 < tm_max and z4 > tm_min) {
                            tm_list[ii] = nr;
                            ii++;
                        }
                    }
                }
            }
        }
        int n_tm_res = ii;

        // 1. Initialization
        for(int nr: range(n_residue)) {
            type_count[nr] = 0;
        }

        int1d surf_cb(n_tm_res, 0);
        for(int ii=0; ii<n_tm_res; ii++){
            int nr = tm_list[ii];
            for (int j=0;j<3;j++) {
                bbx[ii*3+j] = all_pos(0, nr*3+j);
                bby[ii*3+j] = all_pos(1, nr*3+j);
                bbz[ii*3+j] = all_pos(2, nr*3+j);
            }
            cbx[ii] = all_pos(0, index_cb[nr]);
            cby[ii] = all_pos(1, index_cb[nr]);
            cbz[ii] = all_pos(2, index_cb[nr]);
            sfx[nr] = 0.0;
            sfy[nr] = 0.0;
        }

        // 2, calculate the center of XY for TM residues
        float cx = 0.;
        float cy = 0.;
        for(int nr: range(n_tm_res)) {
            cx += cbx[nr];
            cy += cby[nr];
        }
        cx /= n_tm_res;
        cy /= n_tm_res;

        // Start rotating to different angles to identify the surface and force direction
        for(int ai=0; ai<n_rotation; ai++) {

            // 3, Rotate the coordinates
            float angle = ai * dangle;
            RotateXY(cbx, cby, cx, cy, angle, n_tm_res); // for cb atoms of TM residues
            RotateXY(bbx, bby, cx, cy, angle, n_tm_res); // for backbone atoms of TM residues
      
            // 4. Determine boundaries on XY plane
            float minx = 0.0;
            float maxx = 0.0;
            FindMinMax(bbx, minx, maxx, false, n_tm_res*3);
            FindMinMax(cbx, minx, maxx, true,  n_tm_res);
            float miny = 0.0;
            float maxy = 0.0;
            FindMinMax(bby, miny, maxy, false, n_tm_res*3);
            FindMinMax(cby, miny, maxy, true,  n_tm_res);

            minx -= atom_radius;
            maxx += atom_radius;
            miny -= atom_radius;
            maxy += atom_radius;

            // 5. Build the voxels
            float lenx     = maxx-minx;
            float leny     = maxy-miny;
            int   ngrid_x  = (int)((lenx+0.1)/GRID_X)+3;
            int   ngrid_y  = (int)((leny+0.1)/GRID_Y)+3;
            int   ngrid_z  = NGRIDZ;

            float box_lenx = ngrid_x*GRID_X;
            float box_leny = ngrid_y*GRID_Y;

            float3 side    = make_vec3( minx-0.5f*(box_lenx-lenx), miny-0.5f*(box_leny-leny), tm_min-atom_radius );
            float3 bscale1 = make_vec3( ngrid_x/box_lenx, ngrid_y/box_leny, r_GRID_Z );
            float3 bscale2 = make_vec3( SNGRID_X*ngrid_x/box_lenx, SNGRID_Y*ngrid_y/box_leny, SNGRID_Z*r_GRID_Z );
            int3d voxels( ngrid_x, int2d(ngrid_y, int1d(ngrid_z, -1)) );

            // 6. Put backbone and Cb atoms in voxels. Label the voxel as 1
            vector<int3> vids(n_tm_res*4);
            for(int nr: range(n_tm_res)) {
                for(int i=0;i<3;i++) {
                    float3 aa_xyz = make_vec3(bbx[nr*3+i], bby[nr*3+i], bbz[nr*3+i]);
                    int iv = nr*4+i;
                    vids[iv] = vec_floor((aa_xyz - side) * bscale1);
                    voxels[vids[iv].x()][vids[iv].y()][vids[iv].z()] = 1;

                    for(int j=0;j<n_aux;j++) {
                        float dx = 0.0;
                        float dy = auxiliary_points[j].xz[0]*atom_radius;
                        float dz = auxiliary_points[j].xz[1]*atom_radius;
                        float3 aa_aux_xyz = make_vec3(bbx[nr*3+i]+dx, bby[nr*3+i]+dy , bbz[nr*3+i]+dz);
                        int3 vid_aux = vec_floor((aa_aux_xyz - side) * bscale1);
                        voxels[vid_aux.x()][vid_aux.y()][vid_aux.z()] = 1;
                    }
                    for(int j=0;j<n_aux;j++) {
                        float dx = auxiliary_points[j].xz[0]*atom_radius;
                        float dy = 0.0;
                        float dz = auxiliary_points[j].xz[1]*atom_radius;
                        float3 aa_aux_xyz = make_vec3(bbx[nr*3+i]+dx, bby[nr*3+i]+dy , bbz[nr*3+i]+dz);
                        int3 vid_aux = vec_floor((aa_aux_xyz - side) * bscale1);
                        voxels[vid_aux.x()][vid_aux.y()][vid_aux.z()] = 1;
                    }

                }

                float3 cb_xyz = make_vec3(cbx[nr], cby[nr], cbz[nr]);
                int iv = nr*4+3;
                vids[iv] = vec_floor((cb_xyz - side) * bscale1);
                voxels[vids[iv].x()][vids[iv].y()][vids[iv].z()] = 1;

                for(int j=0;j<n_aux;j++) {
                    float dx = 0.0;
                    float dy = auxiliary_points[j].xz[0]*atom_radius;
                    float dz = auxiliary_points[j].xz[1]*atom_radius;
                    float3 cb_aux_xz = make_vec3(cbx[nr]+dx, cby[nr]+dy , cbz[nr]+dz);
                    int3 cb_vid_aux = vec_floor((cb_aux_xz - side) * bscale1);
                    voxels[cb_vid_aux.x()][cb_vid_aux.y()][cb_vid_aux.z()] = 1;
                }
                for(int j=0;j<n_aux;j++) {
                    float dx = auxiliary_points[j].xz[0]*atom_radius;
                    float dy = 0.0;
                    float dz = auxiliary_points[j].xz[1]*atom_radius;
                    float3 cb_aux_yz = make_vec3(cbx[nr]+dx, cby[nr]+dy , cbz[nr]+dz);
                    int3 cb_vid_aux = vec_floor((cb_aux_yz - side) * bscale1);
                    voxels[cb_vid_aux.x()][cb_vid_aux.y()][cb_vid_aux.z()] = 1;
                }
            }

            // 7. find the surface voxels 
            int ngrid_xz = ngrid_x*ngrid_z;
            int ngrid_yz = ngrid_y*ngrid_z;
            int big_int = ngrid_x*ngrid_y*ngrid_z;

            if (method_type == 0) {
                int1d xmin_surf(ngrid_yz, big_int);
                int1d xmax_surf(ngrid_yz, -1);
                int1d ymin_surf(ngrid_xz, big_int);
                int1d ymax_surf(ngrid_xz, -1);
                for(int i=0;i<ngrid_x;i++) {
                    for(int j=0;j<ngrid_y;j++) {
                        for(int k=0;k<ngrid_z;k++) {
                            if (voxels[i][j][k] == 1) {
                                int ik = i*ngrid_z+k;
                                int jk = j*ngrid_z+k;
                                if (i < xmin_surf[jk]) xmin_surf[jk] = i;
                                if (i > xmax_surf[jk]) xmax_surf[jk] = i;
                                if (j < ymin_surf[ik]) ymin_surf[ik] = j;
                                if (j > ymax_surf[ik]) ymax_surf[ik] = j;
                            }
                        }
                    }
                }
                // label the x surface voxels as 10 (left) or 100 (right)
                for(int jk=0;jk<ngrid_yz;jk++) {
                    int j = jk/ngrid_z;
                    int k = jk%ngrid_z;
                    int i = xmin_surf[jk];
                    if (i < big_int) voxels[i][j][k] += 10;
                    i = xmax_surf[jk];
                    if (i > -1) voxels[i][j][k] += 100;
                }
                // label the y surface voxels as 1000 (front) or 10000 (behind)
                for(int ik=0;ik<ngrid_xz;ik++) {
                    int i = ik/ngrid_z;
                    int k = ik%ngrid_z;
                    int j = ymin_surf[ik];
                    if (j < big_int) voxels[i][j][k] += 1000;
                    j = ymax_surf[ik];
                    if (j > -1) voxels[i][j][k] += 10000;
                }
            }
            else if (method_type == 1) {

                // 7.1 Initialize the lipid voxel
                for(int k=0;k<ngrid_z;k++) {
                    for(int i=0;i<ngrid_x;i++) {
                        voxels[i][0][k]         = 0;
                        voxels[i][ngrid_y-1][k] = 0;
                    }

                    for(int j=0;j<ngrid_y;j++) {
                        voxels[0][j][k]         = 0;
                        voxels[ngrid_x-1][j][k] = 0;
                    }
                }
                // 7.2 find the head gropu voxels
                for(int k=0; k<ngrid_z; k += ngrid_z-1) {
                    bool find_all = false;
                    while (not find_all) {
                        int new_find = 0;
                        for(int i=1;i<ngrid_x-1;i++) {
                            for(int j=1;j<ngrid_y-1;j++) {
                                if (voxels[i][j][k] == -1) {
                                    for(int ii=i-1;ii<i+2;ii++) {
                                        for(int jj=j-1;jj<j+2;jj++) {
                                            if (voxels[ii][jj][k] == 0) {
                                                voxels[i][j][k] = 0;
                                                new_find++;       
                                                break;
                                            }
                                        }
                                    }
                                }
                            }
                        }

                        if(new_find==0) find_all = true;
                    }
                }
                // 7.3 find the lower tail group voxels
                for(int k=1; k<ngrid_z/2; k++ ) {
                    bool find_all = false;
                    while (not find_all) {
                        int new_find = 0;
                        for(int i=1;i<ngrid_x-1;i++) {
                            for(int j=1;j<ngrid_y-1;j++) {
                                if (voxels[i][j][k] == -1) {
                                    for(int ii=i-1;ii<i+2;ii++) {
                                        for(int jj=j-1;jj<j+2;jj++) {
                                            if (voxels[ii][jj][k] == 0 and voxels[ii][jj][k-1] == 0) {
                                                voxels[i][j][k] = 0;
                                                new_find++;       
                                                break;
                                            }
                                        }
                                    }
                                }
                            }
                        }

                        if(new_find==0) find_all = true;
                    }
                }
                // 7.4 find the upper tail group voxels
                for(int k=ngrid_z-1; k>=ngrid_z/2; k-- ) {
                    bool find_all = false;
                    while (not find_all) {
                        int new_find = 0;
                        for(int i=1;i<ngrid_x-1;i++) {
                            for(int j=1;j<ngrid_y-1;j++) {
                                if (voxels[i][j][k] == -1) {
                                    for(int ii=i-1;ii<i+2;ii++) {
                                        for(int jj=j-1;jj<j+2;jj++) {
                                            if (voxels[ii][jj][k] == 0 and voxels[ii][jj][k+1] == 0) {
                                                voxels[i][j][k] = 0;
                                                new_find++;       
                                                break;
                                            }
                                        }
                                    }
                                }
                            }
                        }

                        if(new_find==0) find_all = true;
                    }
                }

                // 7.5 find the surface voxels 
                for(int i=1;i<ngrid_x-1;i++) {
                    for(int j=1;j<ngrid_y-1;j++) {
                        for(int k=0;k<ngrid_z;k++) {
                            if (voxels[i][j][k] > 0) {
                                // label the x surface voxels as 10 (left) or 100 (right)
                                // label the y surface voxels as 1000 (left) or 10000 (right)
                                if (voxels[i-1][j][k] == 0) voxels[i][j][k] += 10;
                                if (voxels[i+1][j][k] == 0) voxels[i][j][k] += 100;
                                if (voxels[i][j-1][k] == 0) voxels[i][j][k] += 1000;
                                if (voxels[i][j+1][k] == 0) voxels[i][j][k] += 10000;
                            }
                        }
                    }
                }
            }

            big_int = ngrid_x*SNGRID_X+ngrid_y*SNGRID_Y+ngrid_z*SNGRID_Z;
            // 8. more refined stage: find the interface CB atom in the interface voxel
            int1d xmin_surf_fine(SNGRID_Y*SNGRID_Z*ngrid_yz, big_int);
            int1d xmax_surf_fine(SNGRID_Y*SNGRID_Z*ngrid_yz, -1);
            int1d ymin_surf_fine(SNGRID_X*SNGRID_Z*ngrid_xz, big_int);
            int1d ymax_surf_fine(SNGRID_X*SNGRID_Z*ngrid_xz, -1);
            // search the interface backbone atoms
            for(int nr: range(n_tm_res)) {
                for(int i=0;i<3;i++) {
                    int iv        = nr*4+i;
                    int vtype     = voxels[vids[iv].x()][vids[iv].y()][vids[iv].z()];

                    if (vtype <= 1) continue;
                    float3 aa_xyz = make_vec3(bbx[nr*3+i], bby[nr*3+i], bbz[nr*3+i]);
                    int3 vid2     = vec_floor((aa_xyz - side) * bscale2);

                    int kj        = vid2.y()*ngrid_z*SNGRID_Z+vid2.z();
                    int ki        = vid2.x()*ngrid_z*SNGRID_Z+vid2.z();

                    int v5 = vtype/10000;
                    int v4 = (vtype%10000)/1000;
                    int v3 = ((vtype%10000)%1000)/100;
                    int v2 = (((vtype%10000)%1000)%100)/10;

                    if (v2 == 1 and vid2.x() < xmin_surf_fine[kj]) xmin_surf_fine[kj] = vid2.x();
                    if (v3 == 1 and vid2.x() > xmax_surf_fine[kj]) xmax_surf_fine[kj] = vid2.x();
                    if (v4 == 1 and vid2.y() < ymin_surf_fine[ki]) ymin_surf_fine[ki] = vid2.y();
                    if (v5 == 1 and vid2.y() > ymax_surf_fine[ki]) ymax_surf_fine[ki] = vid2.y();

                    for(int j=0;j<n_aux;j++) {
                        float dx = 0.0;
                        float dy = auxiliary_points[j].xz[0]*atom_radius;
                        float dz = auxiliary_points[j].xz[1]*atom_radius;
                        float3 aa_aux_xyz = make_vec3(bbx[nr*3+i]+dx, bby[nr*3+i]+dy , bbz[nr*3+i]+dz);
                        int3 vid_aux = vec_floor((aa_aux_xyz - side) * bscale2);
                        int kj = vid_aux.y()*ngrid_z*SNGRID_Z+vid2.z();
                        if (v2 == 1 and vid_aux.x() < xmin_surf_fine[kj]) xmin_surf_fine[kj] = vid_aux.x();
                        if (v3 == 1 and vid_aux.x() > xmax_surf_fine[kj]) xmax_surf_fine[kj] = vid_aux.x();
                    }

                    for(int j=0;j<n_aux;j++) {
                        float dx = auxiliary_points[j].xz[0]*atom_radius;
                        float dy = 0.0;
                        float dz = auxiliary_points[j].xz[1]*atom_radius;
                        float3 aa_aux_xyz = make_vec3(bbx[nr*3+i]+dx, bby[nr*3+i]+dy , bbz[nr*3+i]+dz);
                        int3 vid_aux = vec_floor((aa_aux_xyz - side) * bscale2);
                        int ki = vid_aux.x()*ngrid_z*SNGRID_Z+vid2.z();
                        if (v4 == 1 and vid_aux.y() < ymin_surf_fine[ki]) ymin_surf_fine[ki] = vid_aux.y();
                        if (v5 == 1 and vid_aux.y() > ymax_surf_fine[ki]) ymax_surf_fine[ki] = vid_aux.y();
                    }
                }
            }

            // search the interface CB aux atoms
            for(int nr: range(n_tm_res)) {
                int iv        = nr*4+3;
                int vtype     = voxels[vids[iv].x()][vids[iv].y()][vids[iv].z()];

                if (vtype <= 1) continue;
                float3 cb_xyz = make_vec3(cbx[nr], cby[nr], cbz[nr]);
                int3 vid2     = vec_floor((cb_xyz - side) * bscale2);
                int v5        = vtype/10000;
                int v4        = (vtype%10000)/1000;
                int v3        = ((vtype%10000)%1000)/100;
                int v2        = (((vtype%10000)%1000)%100)/10;

                for(int j=0;j<n_aux;j++) {
                    float dx = 0.0;
                    float dy = auxiliary_points[j].xz[0]*atom_radius;
                    float dz = auxiliary_points[j].xz[1]*atom_radius;
                    float3 aa_aux_xyz = make_vec3(cbx[nr]+dx, cby[nr]+dy , cbz[nr]+dz);
                    int3 vid_aux = vec_floor((aa_aux_xyz - side) * bscale2);
                    int kj = vid_aux.y()*ngrid_z*SNGRID_Z+vid2.z();
                    if (v2 == 1 and vid_aux.x() < xmin_surf_fine[kj]) xmin_surf_fine[kj] = vid_aux.x();
                    if (v3 == 1 and vid_aux.x() > xmax_surf_fine[kj]) xmax_surf_fine[kj] = vid_aux.x();
                }

                for(int j=0;j<n_aux;j++) {
                    float dx = auxiliary_points[j].xz[0]*atom_radius;
                    float dy = 0.0;
                    float dz = auxiliary_points[j].xz[1]*atom_radius;
                    float3 aa_aux_xyz = make_vec3(cbx[nr]+dx, cby[nr]+dy , cbz[nr]+dz);
                    int3 vid_aux = vec_floor((aa_aux_xyz - side) * bscale2);
                    int ki = vid_aux.x()*ngrid_z*SNGRID_Z+vid2.z();
                    if (v4 == 1 and vid_aux.y() < ymin_surf_fine[ki]) ymin_surf_fine[ki] = vid_aux.y();
                    if (v5 == 1 and vid_aux.y() > ymax_surf_fine[ki]) ymax_surf_fine[ki] = vid_aux.y();
                }
            }

            // search the interface CB atoms
            int1d xmin_surf_id(SNGRID_Y*SNGRID_Z*ngrid_yz, -1);
            int1d xmax_surf_id(SNGRID_Y*SNGRID_Z*ngrid_yz, -1);
            int1d ymin_surf_id(SNGRID_X*SNGRID_Z*ngrid_xz, -1);
            int1d ymax_surf_id(SNGRID_X*SNGRID_Z*ngrid_xz, -1);
            for(int nr: range(n_tm_res)) {
                int iv        = nr*4+3;
                int vtype     = voxels[vids[iv].x()][vids[iv].y()][vids[iv].z()];

                if (vtype <= 1) continue;
                float3 cb_xyz = make_vec3(cbx[nr], cby[nr], cbz[nr]);
                int3 vid2     = vec_floor((cb_xyz - side) * bscale2);
                int kj        = vid2.y()*ngrid_z*SNGRID_Z+vid2.z();
                int ki        = vid2.x()*ngrid_z*SNGRID_Z+vid2.z();

                int v5 = vtype/10000;
                int v4 = (vtype%10000)/1000;
                int v3 = ((vtype%10000)%1000)/100;
                int v2 = (((vtype%10000)%1000)%100)/10;
                if ( v2 == 1 and vid2.x() <= xmin_surf_fine[kj]) {
                    xmin_surf_fine[kj] = vid2.x();
                    xmin_surf_id[kj]   = nr;
                }
                if ( v3 == 1 and vid2.x() >= xmax_surf_fine[kj]) {
                    xmax_surf_fine[kj] = vid2.x();
                    xmax_surf_id[kj]   = nr;
                }
                if ( v4 == 1 and vid2.y() <= ymin_surf_fine[ki]) {
                    ymin_surf_fine[ki] = vid2.y();
                    ymin_surf_id[ki]   = nr;
                }
                if ( v5 == 1 and vid2.y() >= ymax_surf_fine[ki]) {
                    ymax_surf_fine[ki] = vid2.y();
                    ymax_surf_id[ki]   = nr;
                }
            }

            // count the surface area
            int num_xmin = 0; for (int i:xmin_surf_fine) if (i<big_int) num_xmin++;
            int num_xmax = 0; for (int i:xmax_surf_fine) if (i > -1   ) num_xmax++;
            int num_ymin = 0; for (int i:ymin_surf_fine) if (i<big_int) num_ymin++;
            int num_ymax = 0; for (int i:ymax_surf_fine) if (i > -1   ) num_ymax++;
            // count the number of surface CB atoms
            int num_xmin_cb = 0; for (int i:xmin_surf_id) if (i > -1) num_xmin_cb++;
            int num_xmax_cb = 0; for (int i:xmax_surf_id) if (i > -1) num_xmax_cb++;
            int num_ymin_cb = 0; for (int i:ymin_surf_id) if (i > -1) num_ymin_cb++;
            int num_ymax_cb = 0; for (int i:ymax_surf_id) if (i > -1) num_ymax_cb++;

            float area_xz = GRID_X/SNGRID_X * GRID_Z/SNGRID_Z;
            float area_yz = GRID_Y/SNGRID_Y * GRID_Z/SNGRID_Z;
            float area_xmin = num_xmin*area_yz;
            float area_xmax = num_xmax*area_yz;
            float area_ymin = num_ymin*area_xz;
            float area_ymax = num_ymax*area_xz;

            float area_x = (area_xmin > area_xmax) ? area_xmin : area_xmax;
            float area_y = (area_ymin > area_ymax) ? area_ymin : area_ymax;
            float total_force_x = area_x*FORCE_UNIT;
            float total_force_y = area_y*FORCE_UNIT;

            float fxmin_per_CB = total_force_x/num_xmin_cb;
            float fxmax_per_CB = total_force_x/num_xmax_cb;
            float fymin_per_CB = total_force_y/num_ymin_cb;
            float fymax_per_CB = total_force_y/num_ymax_cb;

            // 9. Record all surface CB atoms. Assign forces derection according to the interface
            for(int nr: range(n_tm_res)) {
                fx[nr] = 0.0;
                fy[nr] = 0.0;
            }
            for(int nr:xmin_surf_id) {
                if (nr >= 0) {
                    surf_cb[nr] = 1; 
                    fx[nr]  = fxmin_per_CB;
                }
            }
            for(int nr:xmax_surf_id) {
                if (nr >= 0) {
                    surf_cb[nr] = 1; 
                    fx[nr] -= fxmax_per_CB;
                }
            }
            for(int nr:ymin_surf_id) {
                if (nr >= 0) {
                    surf_cb[nr] = 1;
                    fy[nr]  = fymin_per_CB;
                }
            }
            for(int nr:ymax_surf_id) {
                if (nr >= 0) {
                    surf_cb[nr] = 1;
                    fy[nr] -= fymax_per_CB;
                }
            }


            // count the 
            for(int ns: range(n_tm_res)) {
                if (surf_cb[ns] == 0) continue;
                int nr = tm_list[ns];
                type_count[nr] += 1;
            }

            RotateXY( fx,  fy, 0.f, 0.f, -angle, n_tm_res);
            RotateXY(cbx, cby,  cx,  cy, -angle, n_tm_res);
            RotateXY(bbx, bby,  cx,  cy, -angle, n_tm_res);


            for(int ns: range(n_tm_res)) {
                int nr   = tm_list[ns];
                if (surf_cb[ns] == 0) continue;
                sfx[nr] += fx[ns];
                sfy[nr] += fy[ns];
            }
        }

        for(int nr: tm_list) {
            if (type_count[nr] <= 0) continue;
            sfx[nr] /= type_count[nr];
            sfy[nr] /= type_count[nr];
        }

        // 9. Assign mempot_switch  
        for(int nr: range(n_residue)) {
            if (type_count[nr] > 0)
                output(0, nr) = type_count[nr]*r_n_rotation;
            else
                output(0, nr) = 0.0;
            output(1, nr) = sfx[nr];
            output(2, nr) = sfy[nr];
        }
    }

    virtual void propagate_deriv() override {
       //pos2.sens(nd,p.index_pos2) += sens(nd+n_dim1,ne);
    }

};
static RegisterNodeType<Surface, 1> surface_node("surface");

