/*
 *  SPDX-FileCopyrightText: 2024 Faris Serdar Taşel <fst@cankaya.edu.tr>
 *  SPDX-FileCopyrightText: 2024 Efe Çiftci <efeciftci@cankaya.edu.tr>
 *
 *  SPDX-License-Identifier: GPL-3.0-or-later
 */

#include "imod.h"

int imod_endian_reverse(int x) {
    return ((x & 0xff000000) >> 24) | ((x & 0x00ff0000) >> 8) |
           ((x & 0x0000ff00) << 8) | ((x & 0x000000ff) << 24);
}

float imod_endian_reverse(float x) {
    void *p = &x;
    int y = imod_endian_reverse(*((int *)p));
    p = &y;
    return *((float *)p);
}

void savePolyArrayAsIMOD(vector<poly25d> &p, int x_size, int y_size, int z_size,
                         double xy_scale, double z_scale, int x_shift,
                         int y_shift) {
    int id;
    uint i, j, k;
    imod_model model;
    imod_objt objt;
    imod_cont cont;
    imod_mesh mesh;
    float pt;

    string outputfn;
    outputfn = (boost::format("%s%s%s.mod") % DESTPATH % "imod_" % FNAME).str();
    outputfn = (boost::format(outputfn) % SLICE_START).str();
    cout << "Saving: " << outputfn << endl;
    ofstream f(outputfn, ios::binary);

    // Write header
    id = IMOD_VALUE(IMOD_FILE_ID);
    f.write(reinterpret_cast<char *>(&id), sizeof(int));
    id = IMOD_VALUE(IMOD_VERSION_ID);
    f.write(reinterpret_cast<char *>(&id), sizeof(int));
    // Write model info
    memset(&model, '\0', sizeof(imod_model));
    model.objsize = IMOD_VALUE((int)p.size());
    model.xmax = IMOD_VALUE(x_size);
    model.ymax = IMOD_VALUE(y_size);
    model.zmax = IMOD_VALUE(z_size);
    model.blacklevel = IMOD_VALUE(0);
    model.whitelevel = IMOD_VALUE(255);
    model.xscale = IMOD_VALUE(1.f);
    model.yscale = IMOD_VALUE(1.f);
    model.zscale = IMOD_VALUE(1.f);
    f.write(reinterpret_cast<char *>(&model), sizeof(imod_model));
    // Write objects' data
    for (k = 0; k < p.size(); k++) {
        id = IMOD_VALUE(IMOD_OBJT_ID);
        f.write(reinterpret_cast<char *>(&id), sizeof(int));
        memset(&objt, '\0', sizeof(imod_objt));
        objt.contsize = IMOD_VALUE((int)p[k].t);
        objt.red = IMOD_VALUE(0.f);
        objt.green = IMOD_VALUE(0.f);
        objt.blue = IMOD_VALUE(1.f);
        objt.symsize = 1;
        objt.meshsize = IMOD_VALUE(1);
        sprintf(objt.name, "obj_%d", k + 1);
        f.write(reinterpret_cast<char *>(&objt), sizeof(imod_objt));
        // Write object's contours
        for (i = 0; i < p[k].t; i++) {
            for (polygon &pl : p[k].slice[i]) {
                id = IMOD_VALUE(IMOD_CONT_ID);
                f.write(reinterpret_cast<char *>(&id), sizeof(int));
                memset(&cont, '\0', sizeof(imod_cont));
                cont.psize = IMOD_VALUE((int)pl.size());
                f.write(reinterpret_cast<char *>(&cont), sizeof(imod_cont));
                // Write contour points
                for (j = 0; j < pl.size(); j++) {
                    pt = (float)(pl[j].x * xy_scale + x_shift);
                    pt = IMOD_VALUE(pt);
                    f.write(reinterpret_cast<char *>(&pt), sizeof(float));
                    pt = (float)(y_size - 1 - pl[j].y * xy_scale - y_shift);
                    pt = IMOD_VALUE(pt);
                    f.write(reinterpret_cast<char *>(&pt), sizeof(float));
                    pt = (float)((p[k].start_z + i) * z_scale);
                    pt = IMOD_VALUE(pt);
                    f.write(reinterpret_cast<char *>(&pt), sizeof(float));
                }
            }
        }
        // Get mesh info
        vector<vector<int3>> faces;
        vector<uint> n_faces;
        faces.resize(p[k].t - 1);
        n_faces.resize(p[k].t - 1);
        int n_vert = p[k].slice[0][0].size();
        int n_totalf = 0;

        for (i = 0; i < p[k].t - 1; i++) {
            n_faces[i] =
                getTriangles(p[k].slice[i][0], p[k].slice[i + 1][0], faces[i]);
            n_vert += p[k].slice[i + 1][0].size();
            n_totalf += n_faces[i];
        }
        // Write vertices
        id = IMOD_VALUE(IMOD_MESH_ID);
        f.write(reinterpret_cast<char *>(&id), sizeof(int));
        memset(&mesh, '\0', sizeof(imod_mesh));
        mesh.vsize = IMOD_VALUE(n_vert);
        mesh.lsize = IMOD_VALUE(n_totalf * 5 + 1);
        f.write(reinterpret_cast<char *>(&mesh), sizeof(imod_mesh));
        for (i = 0; i < p[k].t; i++) {
            for (j = 0; j < p[k].slice[i][0].size(); j++) {
                pt = (float)(p[k].slice[i][0][j].x * xy_scale + x_shift);
                pt = IMOD_VALUE(pt);
                f.write(reinterpret_cast<char *>(&pt), sizeof(float));
                pt = (float)(y_size - 1 - p[k].slice[i][0][j].y * xy_scale -
                             y_shift);
                pt = IMOD_VALUE(pt);
                f.write(reinterpret_cast<char *>(&pt), sizeof(float));
                pt = (float)((p[k].start_z + i) * z_scale);
                pt = IMOD_VALUE(pt);
                f.write(reinterpret_cast<char *>(&pt), sizeof(float));
            }
        }
        // Write triangles
        int count = 0;
        for (i = 0; i < p[k].t - 1; i++) {
            for (j = 0; j < n_faces[i]; j++) {
                id = IMOD_VALUE(-25);
                f.write(reinterpret_cast<char *>(&id), sizeof(int));
                id = IMOD_VALUE(faces[i][j].index[0] + count);
                f.write(reinterpret_cast<char *>(&id), sizeof(int));
                id = IMOD_VALUE(faces[i][j].index[1] + count);
                f.write(reinterpret_cast<char *>(&id), sizeof(int));
                id = IMOD_VALUE(faces[i][j].index[2] + count);
                f.write(reinterpret_cast<char *>(&id), sizeof(int));
                id = IMOD_VALUE(-22);
                f.write(reinterpret_cast<char *>(&id), sizeof(int));
            }
            count += p[k].slice[i][0].size();
        }
        id = IMOD_VALUE(-1);
        f.write(reinterpret_cast<char *>(&id), sizeof(int));
    }
    // EOF
    id = IMOD_VALUE(IMOD_IEOF_ID);
    f.write(reinterpret_cast<char *>(&id), sizeof(int));
    f.close();
}
