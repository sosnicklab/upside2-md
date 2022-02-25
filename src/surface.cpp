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
#define GRID_X 8
#define GRID_Y 8
#define GRID_Z 6
#define SNGRID_X 2
#define SNGRID_Y 2
#define SNGRID_Z 2
#define ForceNorm 0.70710678
#define FORCE_UNIT 2.428602e-5

using namespace h5;
using namespace std;

typedef vector<int>   int1d;
typedef vector<int1d> int2d;
typedef vector<int2d> int3d;
typedef map<int, int1d> map_int1d;

void RotateXY(vector<float>& x, vector<float>& y, float rx, float ry, float angle) {
    int n = x.size();
    for (int na=0; na< n; ++na ) {
        float dx = x[na] - rx;
        float dy = y[na] - ry;
        if (dx==0 and dy==0) continue;
        x[na]    = rx + cos(angle) * dx - sin(angle) * dy;
        y[na]    = ry + sin(angle) * dx + cos(angle) * dy;
    }
}

void FindMinMax(vector<float>& X, float& minx, float& maxx, bool initialization ) {
    for (float& x: X) {
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
    }
}

struct Surface : public CoordNode
{
    int n_residue;  // number of residues
    int n_used_res;
    int n_rotation; // The number of rotations, from 0 to 90 degrees
    float r_n_rotation; // The number of rotations, from 0 to 90 degrees
    float dangle;   // the rotation angle of every rotation

    CoordNode& pos;  // N, CA, CB C atoms

    int1d index_cb;
    int1d type_count;
    int1d used_list;

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

    //read_attribute<int>(grp, ".", "n_dim1")+read_attribute<int>(grp, ".", "n_dim2")),
    Surface(hid_t grp, CoordNode& pos_):
        CoordNode( get_dset_size(1,grp,"index_cb")[0], 3 ),
        n_residue( get_dset_size(1,grp,"index_cb")[0] ),
        n_used_res(n_residue),
        n_rotation( read_attribute<int>(grp, ".", "n_rotation") ),
        r_n_rotation( 1./n_rotation ),
        dangle(0.5*M_PI_F*r_n_rotation),

        pos(pos_),
        index_cb(n_residue),
        type_count(n_residue, -1),
        used_list(n_residue),

        cbx(n_residue),
        cby(n_residue),
        cbz(n_residue),
        bbx(n_residue*3),
        bby(n_residue*3),
        bbz(n_residue*3),
        fx(n_residue,  0.),
        fy(n_residue,  0.),
        sfx(n_residue, 0.),
        sfy(n_residue, 0.)
    {
        check_elem_width_lower_bound(pos, 3);
        traverse_dset<1,  int>(grp, "index_cb", [&](size_t rt, index_t x){ index_cb[rt] = x;});

        // 0, no exclusion; 
        // 1, exclusion by exclusion_list; 
        // 2, exclusion by z
        int exclude_way = read_attribute<int>(grp, ".", "exclusion"); 

        if (exclude_way == 0) {
            for(int nr: range(n_residue)) used_list[nr] = nr;
        }
        else if  (exclude_way == 1) {
            n_used_res = get_dset_size(1,grp,"included_list")[0];
            used_list.resize(n_used_res);
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

        int1d surf_cb(n_used_res, 0);

        // 0. Initialization
        int ii = 0;
        for(int nr: used_list) {
            for (int j=0;j<3;j++) {
                bbx[ii*3+j] = all_pos(0, nr*3+j);
                bby[ii*3+j] = all_pos(1, nr*3+j);
                bbz[ii*3+j] = all_pos(2, nr*3+j);
            }
            cbx[ii] = all_pos(0, index_cb[nr]);
            cby[ii] = all_pos(1, index_cb[nr]);
            cbz[ii] = all_pos(2, index_cb[nr]);
            ii++;

            sfx[nr] = 0.0;
            sfy[nr] = 0.0;
            type_count[nr] = 0;
        }

        // 1, calculate the center of XY for TM residues
        float cx = 0.;
        float cy = 0.;
        for(int nr: range(n_used_res)) {
            cx += cbx[nr];
            cy += cby[nr];
        }
        cx /= n_used_res;
        cy /= n_used_res;

        // 2, find the minz and maxz
        float minz = 0.f;
        float maxz = 0.f;
        FindMinMax(bbz, minz, maxz, false);
        FindMinMax(cbz, minz, maxz, true);

        // Start rotating to different angles to identify the surface and force direction
        for(int ai=0; ai<n_rotation; ai++) {

            // 3, Rotate the coordinates
            float angle = ai * dangle;
            RotateXY(cbx, cby, cx, cy, angle); // for cb atoms of TM residues
            RotateXY(bbx, bby, cx, cy, angle); // for backbone atoms of TM residues
      
            // 4. Determine boundaries on XY plane
            float minx = 0.0;
            float maxx = 0.0;
            FindMinMax(bbx, minx, maxx, false);
            FindMinMax(cbx, minx, maxx, true);
            float miny = 0.0;
            float maxy = 0.0;
            FindMinMax(bby, miny, maxy, false);
            FindMinMax(cby, miny, maxy, true);

            // 5. Build the voxels
            float lenx     = maxx-minx;
            float leny     = maxy-miny;
            float lenz     = maxz-minz;
            int   ngrid_x  = (int)((lenx+0.1)/GRID_X)+1;
            int   ngrid_y  = (int)((leny+0.1)/GRID_Y)+1;
            int   ngrid_z  = (int)((lenz+0.1)/GRID_Z)+1;
            float box_lenx = ngrid_x*GRID_X;
            float box_leny = ngrid_y*GRID_Y;
            float box_lenz = ngrid_z*GRID_Z;
            float3 side    = make_vec3( minx-0.5f*(box_lenx-lenx), miny-0.5f*(box_leny-leny), minz-0.5f*(box_lenz-lenz) );
            float3 bscale1 = make_vec3( ngrid_x/box_lenx, ngrid_y/box_leny, ngrid_z/box_lenz );
            float3 bscale2 = make_vec3( SNGRID_X*ngrid_x/box_lenx, SNGRID_Y*ngrid_y/box_leny, SNGRID_Z*ngrid_z/box_lenz );
            int3d voxels( ngrid_x, int2d(ngrid_y, int1d(ngrid_z,0)) );

            // 6. Put backbone and Cb atoms in voxels. Label the voxel as 1
            vector<int3> vids(n_used_res*4);
	    for(int nr: range(n_used_res)) {
                for(int i=0;i<3;i++) {
                    float3 aa_xyz = make_vec3(bbx[nr*3+i], bby[nr*3+i], bbz[nr*3+i]);
                    int iv = nr*4+i;
                    vids[iv] = vec_floor((aa_xyz - side) * bscale1);
                    voxels[vids[iv].x()][vids[iv].y()][vids[iv].z()] = 1;
                    //int3 vid      = vec_floor((aa_xyz - side) * bscale1);
                    //voxels[vid.x()][vid.y()][vid.z()] = 1;
                }
                float3 cb_xyz = make_vec3(cbx[nr], cby[nr], cbz[nr]);
                //int3 vid      = vec_floor((cb_xyz - side) * bscale1);
                //voxels[vid.x()][vid.y()][vid.z()] = 1;
                int iv = nr*4+3;
                vids[iv] = vec_floor((cb_xyz - side) * bscale1);
                voxels[vids[iv].x()][vids[iv].y()][vids[iv].z()] = 1;
            }

            // 7. find the surface voxels 
            int ngrid_xz = ngrid_x*ngrid_z;
            int ngrid_yz = ngrid_y*ngrid_z;
            int big_int = ngrid_x*ngrid_y*ngrid_z;
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
            // label the x surface voxels as 2 (left) or 3 (right)
            for(int jk=0;jk<ngrid_yz;jk++) {
                int j = jk/ngrid_z;
                int k = jk%ngrid_z;
                int i = xmin_surf[jk];
                if (i < big_int) voxels[i][j][k] += 10;
                i = xmax_surf[jk];
                if (i > -1) voxels[i][j][k] += 100;
            }
            // label the y surface voxels as 4 (front) or  (behind)
            for(int ik=0;ik<ngrid_xz;ik++) {
                int i = ik/ngrid_z;
                int k = ik%ngrid_z;
                int j = ymin_surf[ik];
                if (j < big_int) voxels[i][j][k] += 1000;
                j = ymax_surf[ik];
                if (j > -1) voxels[i][j][k] += 10000;
            }

            big_int = ngrid_x*SNGRID_X+ngrid_y*SNGRID_Y+ngrid_z*SNGRID_Z;
            // 8. more refined stage: find the interface CB atom in the interface voxel
            int1d xmin_surf_fine(SNGRID_Y*SNGRID_Z*ngrid_yz, big_int);
            int1d xmax_surf_fine(SNGRID_Y*SNGRID_Z*ngrid_yz, -1);
            int1d ymin_surf_fine(SNGRID_X*SNGRID_Z*ngrid_xz, big_int);
            int1d ymax_surf_fine(SNGRID_X*SNGRID_Z*ngrid_xz, -1);
            // search the interface backbone atoms
            for(int nr: range(n_used_res)) {
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
                }
            }
            // search the interface CB atoms
            int1d xmin_surf_id(SNGRID_Y*SNGRID_Z*ngrid_yz, -1);
            int1d xmax_surf_id(SNGRID_Y*SNGRID_Z*ngrid_yz, -1);
            int1d ymin_surf_id(SNGRID_X*SNGRID_Z*ngrid_xz, -1);
            int1d ymax_surf_id(SNGRID_X*SNGRID_Z*ngrid_xz, -1);
            for(int nr: range(n_used_res)) {
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
            for(int nr: range(n_used_res)) {
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
            for(int ns: range(n_used_res)) {
                if (surf_cb[ns] == 0) continue;
                int nr = used_list[ns];
                type_count[nr] += 1;
            }

            RotateXY( fx,  fy, 0.f, 0.f, -angle);
            RotateXY(cbx, cby,  cx,  cy, -angle);
            RotateXY(bbx, bby,  cx,  cy, -angle);

            for(int ns: range(n_used_res)) {
                int nr   = used_list[ns];
                if (surf_cb[ns] == 0) continue;
                sfx[nr] += fx[ns];
                sfy[nr] += fy[ns];
            }
        }

        for(int nr: used_list) {
            if (type_count[nr] <= 0) continue;
            sfx[nr] /= type_count[nr];
            sfy[nr] /= type_count[nr];
        }

        // 9. Assign mempot_switch  
        for(int nr: range(n_residue)) {
            if (type_count[nr] > 0)
                output(0, nr) = type_count[nr]*r_n_rotation;
            else
                output(0, nr) = type_count[nr];
            output(1, nr) = sfx[nr];
            output(2, nr) = sfy[nr];
        }
    }

    virtual void propagate_deriv() override {
       //pos2.sens(nd,p.index_pos2) += sens(nd+n_dim1,ne);
    }

};
static RegisterNodeType<Surface, 1> surface_node("surface");

