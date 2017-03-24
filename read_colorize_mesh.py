# coding=utf-8
__author__ = 'Anastasia Bazhutina'
colorize_mesh = 'color(2).ply'
file_result_name = 'color_points(2).vtk'
import numpy as np

def write_points_in_vtk(points, vectors):
    """Запись в файл координат точек, скаляров для каждой точки, векторов в каждой точке

     Parameters
     ----------
     points :
         points - массив точек
    vectors :
        vectors - массив векторов
    """

    len1 = len(points)
    len2 = len(vectors)
    len3 = 0

    print 'start write to file'
    with open(file_result_name, 'w') as out_mesh:
        out_mesh.write("# vtk DataFile Version 2.0\n")
        out_mesh.write("Really cool data\n")
        out_mesh.write("ASCII\n")
        out_mesh.write("DATASET UNSTRUCTURED_GRID\n")
        out_mesh.write("POINTS {pcount} double\n".format(pcount=len1))
        for i in xrange(0, len1):
            out_mesh.write("{0} {1} {2}\n".format(points[i][0], points[i][1], points[i][2]))

        out_mesh.write("CELLS {pcount} {pcount2}\n".format(pcount=len1, pcount2=(len1 * 2)))
        for j in xrange(len1):
            out_mesh.write("1 {0}\n".format(j))
        out_mesh.write("CELL_TYPES {pcount}\n".format(pcount=len1))

        for k in xrange(len1):
            out_mesh.write("1\n")


        if len2 != 0:
            if len3 == 0:
                out_mesh.write("POINT_DATA {pcount}\n".format(pcount=len2))
            out_mesh.write("VECTORS DTMRIFiberOrientation double\n")
            for i in xrange(len2):
                out_mesh.write("{0} {1} {2}\n".format(vectors[i][0], vectors[i][1], vectors[i][2]))



def main():
    f = open(colorize_mesh, 'r')
    n = 0
    for line in f:
        if n == 3:
            N = int(line.split()[2])
            coords = []

            print N
        n += 1
        if n >= 17 and n <= 16 + N:
            s = line.split()
            #print s
            red = int(s[6])
            green = int(s[7])
            blue = int(s[8])
            if not((0 <= red <= 50) and (0 <= green <= 50) and blue != 0):
                coords.append([float(s[0]), float(s[1]), float(s[2])])

    #print coords
    write_points_in_vtk(coords, [])
if __name__ == "__main__":
    main()
