import pygame as pg
import numpy as np
from numba import njit

SCREEN_W, SCREEN_H = 800, 600
FOV = np.pi/2
FOV_V = FOV*SCREEN_H/SCREEN_W


def main():
    pg.init()
    screen = pg.display.set_mode((SCREEN_W, SCREEN_H))
    running = True
    clock = pg.time.Clock()
    surf = pg.surfarray.make_surface(np.zeros((SCREEN_W,SCREEN_H, 3)))

    #points, triangles =  read_obj('VideoShip.obj')
    points = np.asarray([[1, 1, 1, 1, 1], [4, 2, 0, 1, 1], [1, .5, 3, 1, 1]])
    triangles = np.asarray([[0,1,2], [2,1,0]])

    camera = np.asarray([13, 0.5, 2, 3.3, 0])
    light_direction = np.asarray([0, 1, 1])
    
    z_order, TrianglesToRaster, colors = np.zeros(len(triangles)), triangles.copy(), np.zeros((len(triangles), 3))

    while running:
        pg.mouse.set_pos(SCREEN_W/2, SCREEN_H/2)
        surf.fill([50,127,200])
        elapsed_time = clock.tick()*0.001
        for event in pg.event.get():
            if event.type == pg.QUIT:
                running = False
            if event.type == pg.KEYDOWN:
                if event.key == pg.K_ESCAPE:
                    running = False
        
        points = project_points(points, camera)
        z_order, TrianglesToRaster, colors = sort_triangles(points, triangles, camera, light_direction, z_order, TrianglesToRaster, colors) 
        
        for index in np.argsort(np.asarray(z_order)):
            if z_order[index] == 9999:
                break
            p0 = points[TrianglesToRaster[index][0]][3:]
            p1 = points[TrianglesToRaster[index][1]][3:]
            p2 = points[TrianglesToRaster[index][2]][3:]
            pg.draw.polygon(surf, colors[index], [p0, p1, p2])                

        screen.blit(surf, (0,0))
        pg.display.update()
        pg.display.set_caption(str(round(1/(elapsed_time+1e-16), 1)) + str(camera))
        camera = movement(camera, elapsed_time*10)

    return 0

def movement(camera, elapsed_time):
    pressed_keys = pg.key.get_pressed()
    x, y, z, rot, rotv = camera
    diag =  0
    if pg.mouse.get_focused():
        p_mouse = pg.mouse.get_pos()
        rot = (rot + np.clip((p_mouse[0]-SCREEN_W/2)/SCREEN_W, -0.2, .2))%(2*np.pi)
        rotv = rotv + np.clip((p_mouse[1]-SCREEN_H/2)/SCREEN_H, -0.2, .2)
        rotv = np.clip(rotv, -1, 1)

    if pressed_keys[ord('e')]:
        y += elapsed_time
    elif pressed_keys[ord('q')]:
        y -= elapsed_time
        
    if pressed_keys[pg.K_UP] or pressed_keys[ord('w')]:
        x, z, diag = x + elapsed_time*np.cos(rot), z + elapsed_time*np.sin(rot), 1

    elif pressed_keys[pg.K_DOWN] or pressed_keys[ord('s')]:
        x, z, diag = x - elapsed_time*np.cos(rot), z - elapsed_time*np.sin(rot), 1
        
    if pressed_keys[pg.K_LEFT] or pressed_keys[ord('a')]:
        elapsed_time = elapsed_time/(diag+1)
        x, z = x + elapsed_time*np.sin(rot), z - elapsed_time*np.cos(rot)
        
    elif pressed_keys[pg.K_RIGHT] or pressed_keys[ord('d')]:
        elapsed_time = elapsed_time/(diag+1)
        x, z = x - elapsed_time*np.sin(rot), z + elapsed_time*np.cos(rot)

    return np.asarray([x, y, z, rot, rotv])

@njit()
def project_points(points, camera):

    for point in points:

        # Calculate xy angle of vector from camera point to projection point
        h_angle_camera_point = np.arctan((point[2]-camera[2])/(point[0]-camera[0]))
        
        # Check if it isn't pointing backwards
        if abs(camera[0]+np.cos(h_angle_camera_point)-point[0]) > abs(camera[0]-point[0]):
            h_angle_camera_point = (h_angle_camera_point - np.pi)%(2*np.pi)

        # Calculate difference between camera angle and pointing angle
        h_angle = (h_angle_camera_point-camera[3])%(2*np.pi)
        
        # Bring to -pi to pi range
        if h_angle > np.pi:
            h_angle =  h_angle - 2*np.pi
        
        point[3] = SCREEN_W*h_angle/FOV + SCREEN_W/2

        # Calculate xy distance from camera point to projection point
        h_distance = np.sqrt((point[0]-camera[0])**2 + (point[2]-camera[2])**2)
        sine = max(-1, min(1, (camera[1]-point[1])/h_distance))
        
        # Calculate angle to xy plane
        v_angle_camera_point = np.arcsin(sine)

        # Calculate difference between camera verticam angle and pointing vertical angle
        v_angle = (v_angle_camera_point - camera[4])%(2*np.pi)

        if v_angle > np.pi:
            v_angle =  v_angle - 2*np.pi

        point[4] = SCREEN_H*v_angle/FOV_V + SCREEN_H/2

    return points

@njit()
def dot_3d(arr1, arr2):
    return arr1[0]*arr2[0] + arr1[1]*arr2[1] + arr1[2]*arr2[2]

@njit()
def sort_triangles(points, triangles, camera, light_direction, z_order, TrianglesToRaster, colors):
    for i in range(len(triangles)):
        triangle = triangles[i]
        point00 = points[triangle[0]][:3]
        point11 = points[triangle[1]][:3] 
        point22 = points[triangle[2]][:3]

        CameraRay = point00 - camera[:3]
        CameraRay = CameraRay/np.sqrt(CameraRay[0]*CameraRay[0] + CameraRay[1]*CameraRay[1] + CameraRay[2]*CameraRay[2])

        # Use Cross-Product to get surface normal
        line1 = point11 - point00
        line2 = point22 - point00

        normal = np.cross(line1, line2)
        normal = normal/np.sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2])
        
        # get projected 2d points
        point0 = points[triangle[0]][3:5]
        point1 = points[triangle[1]][3:5] 

        minx = min(point0[0],  point1[0],  point1[0])
        maxx = max(point0[0],  point1[0],  point1[0])
        miny = min(point0[1],  point1[1],  point1[1])
        maxy = max(point0[1],  point1[1],  point1[1])

        # check valid values
        if (dot_3d(normal, CameraRay) < 0 and minx > - 1*SCREEN_W and maxx < 2*SCREEN_W
            and miny > - 1*SCREEN_H and maxy < 2*SCREEN_H):
            
            dp = min(1, max(0.1, dot_3d(light_direction, normal)))
            colors[i] = (np.ones(3)*dp*255)

            midpoint = [
                        np.mean(np.asarray([point00[0], point11[0], point22[0]])),
                        np.mean(np.asarray([point00[1], point11[1], point22[1]])),
                        np.mean(np.asarray([point00[2], point11[2], point22[2]]))
                        ]
            zz = np.asarray(midpoint) - camera[:3]
            z_order[i] = -np.sqrt(zz[0]*zz[0] + zz[1]*zz[1] + zz[2]*zz[2])
            TrianglesToRaster[i] = triangle
        else: z_order[i] = 9999
    return z_order, TrianglesToRaster, colors

def read_obj(fileName):
    vertices = []
    triangles = []
    f = open(fileName)
    for line in f:
        if line[:2] == "v ":
            index1 = line.find(" ") + 1
            index2 = line.find(" ", index1 + 1)
            index3 = line.find(" ", index2 + 1)
            
            vertex = [float(line[index1:index2]), float(line[index2:index3]), float(line[index3:-1]), 1, 1]
            #vertex = [float(line[index1:index2]), float(line[index3:-1]), float(line[index2:index3]), 1, 1]
            vertices.append(vertex)

        elif line[0] == "f":
            index1 = line.find(" ") + 1
            index2 = line.find(" ", index1 + 1)
            index3 = line.find(" ", index2 + 1)

            triangles.append([int(line[index1:index2]) - 1, int(line[index2:index3]) - 1, int(line[index3:-1]) - 1])

    f.close()

    return np.asarray(vertices), np.asarray(triangles)


if __name__ == '__main__':
    main()
    pg.quit()
