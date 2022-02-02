import pygame as pg
import numpy as np
from numba import njit

SCREEN_W, SCREEN_H = 800, 600
FOV_V = np.pi/4 # 45 degrees vertical fov
FOV = FOV_V*SCREEN_W/SCREEN_H


def main():
    pg.init()
    screen = pg.display.set_mode((SCREEN_W, SCREEN_H))
    running = True
    clock = pg.time.Clock()
    surf = pg.surfarray.make_surface(np.zeros((SCREEN_W,SCREEN_H, 3)))

    # points = np.asarray([[1, 1, 1, 1, 1], [4, 2, 0, 1, 1], [1, .5, 3, 1, 1]])
    # triangles = np.asarray([[0,1,2], [2,1,0]])
    points, triangles =  read_obj('VideoShip.obj')

    camera = np.asarray([13, 0.5, 2, 3.3, 0])
    light_direction = np.asarray([0, 1, 1])
    
    z_order, ToRaster, colors = np.zeros(len(triangles)), triangles.copy(), np.zeros((len(triangles), 3))

    while running:
        pg.mouse.set_pos(SCREEN_W/2, SCREEN_H/2)
        surf.fill([50,127,200])
        elapsed_time = clock.tick()*0.001
        for event in pg.event.get():
            if event.type == pg.QUIT: running = False
            if event.type == pg.KEYDOWN and event.key == pg.K_ESCAPE: running = False
        
        points = project_points(points, camera)
        z_order, ToRaster, colors = sort_triangles(points, triangles, camera, light_direction, z_order, ToRaster, colors) 
        
        for index in np.argsort(np.asarray(z_order)):
            if z_order[index] == 9999: break
            triangle = [points[ToRaster[index][0]][3:], points[ToRaster[index][1]][3:], points[ToRaster[index][2]][3:]]
            pg.draw.polygon(surf, colors[index], triangle)                

        screen.blit(surf, (0,0)); pg.display.update()
        pg.display.set_caption(str(round(1/(elapsed_time+1e-16), 1)) + ' ' + str(camera))
        camera = movement(camera, elapsed_time*10)

    return 0

def movement(camera, elapsed_time):
    
    x, y, z, rot, rotv = camera
    diag =  0

    if pg.mouse.get_focused():
        p_mouse = pg.mouse.get_pos()
        rot = (rot + np.clip((p_mouse[0]-SCREEN_W/2)/SCREEN_W, -0.2, .2))%(2*np.pi)
        rotv = rotv + np.clip((p_mouse[1]-SCREEN_H/2)/SCREEN_H, -0.2, .2)
        rotv = np.clip(rotv, -.5, .5)
    
    pressed_keys = pg.key.get_pressed()
    if pressed_keys[ord('e')]: y += elapsed_time
    elif pressed_keys[ord('q')]: y -= elapsed_time
        
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
        if h_angle > np.pi: h_angle =  h_angle - 2*np.pi
        
        point[3] = SCREEN_W*h_angle/FOV + SCREEN_W/2

        # Calculate xy distance from camera point to projection point
        h_distance = np.sqrt((point[0]-camera[0])**2 + (point[2]-camera[2])**2)
        sine = max(-1, min(1, (camera[1]-point[1])/h_distance))
        
        # Calculate angle to xy plane
        v_angle_camera_point = np.arcsin(sine)

        # Calculate difference between camera verticam angle and pointing vertical angle
        v_angle = (v_angle_camera_point - camera[4])%(2*np.pi)

        if v_angle > np.pi: v_angle =  v_angle - 2*np.pi

        point[4] = SCREEN_H*v_angle/FOV_V + SCREEN_H/2

    return points

@njit()
def dot_3d(arr1, arr2):
    return arr1[0]*arr2[0] + arr1[1]*arr2[1] + arr1[2]*arr2[2]

@njit()
def sort_triangles(points, triangles, camera, light_direction, z_order, ToRaster, colors):
    for i in range(len(triangles)):
        triangle = triangles[i]

        # Use Cross-Product to get surface normal
        line1 = points[triangle[1]][:3]  - points[triangle[0]][:3]
        line2 = points[triangle[2]][:3] - points[triangle[0]][:3]

        # backface culling with dot product between normal and camera ray
        normal = np.cross(line1, line2)
        normal = normal/np.sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2])
        
        CameraRay = points[triangle[0]][:3] - camera[:3]
        CameraRay = CameraRay/np.sqrt(CameraRay[0]*CameraRay[0] + CameraRay[1]*CameraRay[1] + CameraRay[2]*CameraRay[2])

        # get projected 2d points for filtering of offscreen triangles
        xxs = np.asarray([points[triangle[0]][3:5][0],  points[triangle[1]][3:5][0],  points[triangle[2]][3:5][0]])
        yys = np.asarray([points[triangle[0]][3:5][1],  points[triangle[1]][3:5][1],  points[triangle[2]][3:5][1]])

        # check valid values
        if (dot_3d(normal, CameraRay) < 0   and np.min(xxs) > - SCREEN_W and np.max(xxs) < 2*SCREEN_W
                                            and np.min(yys) > - SCREEN_H and np.max(yys) < 2*SCREEN_H):
            
            # calculate triangle shading
            dp = min(1, max(0.1, dot_3d(light_direction, normal)))
            colors[i] = (np.ones(3)*dp*255)

            # calculate distance to camera for drawing order, negative for inverse order
            dist2cam = (points[triangle[0]][:3] + points[triangle[1]][:3] + points[triangle[2]][:3])/3 - camera[:3]
            z_order[i] = -np.sqrt(dist2cam[0]*dist2cam[0] + dist2cam[1]*dist2cam[1] + dist2cam[2]*dist2cam[2])

            # save current triangle for rasterization
            ToRaster[i] = triangle
        
        # big value for last positions in sort
        else: z_order[i] = 9999

    return z_order, ToRaster, colors

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
