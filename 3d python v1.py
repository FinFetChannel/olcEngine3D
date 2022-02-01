import pygame as pg
import numpy as np

SCREEN_W, SCREEN_H = 1280, 720
FOV = np.pi/2
FOV_V = FOV*SCREEN_H/SCREEN_W


def main():
    pg.init()
    screen = pg.display.set_mode((SCREEN_W, SCREEN_H, 3))
    running = True
    clock = pg.time.Clock()
    surf = pg.surfarray.make_surface(np.zeros((SCREEN_W,SCREEN_H)))

    points, triangles =  read_obj('VideoShip.obj')
    # points = [[1, 1, 0, 1, 1], [1, 2, 0, 1, 1], [1, .5, 1, 1, 1]]
    # triangles = [[0,1,2]]

    camera = np.asarray([0.5, 0.5, 0.5, 4.86, 0])
    light_direction = np.asarray([0, 1, 1])

    while running:
        surf.fill([50,127,200])
        elapsed_time = clock.tick()*0.001
        for event in pg.event.get():
            if event.type == pg.QUIT:
                running = False
            if event.type == pg.KEYDOWN:
                if event.key == pg.K_ESCAPE:
                    running = False
        
        points = project_points(points, camera)
        # camera_vector = np.asarray([camera[0]+np.cos(camera[3]), camera[1]+np.sin(camera[3]), camera[2]+np.sin(camera[4]) ])
        # camera_vector = camera_vector/np.linalg.norm(camera_vector)

        z_order, TrianglesToRaster, colors = [], [], []
        for triangle in triangles:
            point00 = points[triangle[0]][:3]
            point11 = points[triangle[1]][:3] 
            point22 = points[triangle[2]][:3]

            CameraRay = point00 - camera[:3]
            CameraRay = CameraRay/np.linalg.norm(CameraRay)

            # Use Cross-Product to get surface normal
            line1 = point11 - point00
            line2 = point22 - point00

            normal = np.cross(line1, line2)
            normal = normal/np.linalg.norm(normal)
            
            # get projected 2d points
            point0 = points[triangle[0]][3:5]
            point1 = points[triangle[1]][3:5] 
            point2 = points[triangle[2]][3:5]

            minx = min(point0[0],  point1[0],  point1[0])
            maxx = max(point0[0],  point1[0],  point1[0])
            miny = min(point0[1],  point1[1],  point1[1])
            maxy = max(point0[1],  point1[1],  point1[1])

            # check valid values
            if (np.dot(normal, CameraRay) < 0 and minx > - 0.3*SCREEN_W and maxx < 1.3*SCREEN_W
                and miny > - 0.3*SCREEN_H and maxy < 1.3*SCREEN_H):
                
                dp = min(1, max(0.1, np.dot(light_direction, normal)))
                #print(dp)
                colors.append(np.ones(3)*dp*255)

                midpoint = [
                            np.mean([point00[0], point11[0], point22[0]]),
                            np.mean([point00[1], point11[1], point22[1]]),
                            np.mean([point00[2], point11[2], point22[2]])
                            ]
                z_order.append( -np.linalg.norm(np.asarray(midpoint) - np.asarray(camera[:3])))
                TrianglesToRaster.append([point0, point1, point2])
        
        for index in np.argsort(np.asarray(z_order)):
            pg.draw.polygon(surf, colors[index], TrianglesToRaster[index])                

        screen.blit(surf, (0,0))
        pg.display.update()
        pg.display.set_caption(str(round(1/(elapsed_time+1e-16), 1)))# + str(camera))
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
        rotv = np.clip(rotv, -3, 3)
        pg.mouse.set_pos(SCREEN_W/2, SCREEN_H/2)

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

def project_points(points, camera):

    for point in points:
        # camera_vector = np.asarray([camera[0]+np.cos(camera[3]), camera[1]+np.sin(camera[3]), camera[2]+np.sin(camera[4]) ])
        # camera_vector = camera_vector/np.linalg.norm(camera_vector)

        # point_vector = np.asarray([point[0] - camera[0], point[1] - camera[1], point[2] - camera[2]])
        # point_vector = point_vector/np.linalg.norm(point_vector)

        # if np.dot(camera_vector, point_vector) < -0.5:
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

        # else:
        #     point[3] = -SCREEN_W
        #     point[4] = -SCREEN_H

    return points

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
