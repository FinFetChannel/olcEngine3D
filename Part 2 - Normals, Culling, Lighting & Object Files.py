import pygame as pg
import numpy as np

SCREEN_WIDTH, SCREEN_HEIGHT = 800, 600

pg.init()
screen = pg.display.set_mode((SCREEN_WIDTH,SCREEN_HEIGHT))
running = True
clock = pg.time.Clock()
surf = pg.surfarray.make_surface(np.zeros((SCREEN_WIDTH,SCREEN_HEIGHT, 3)))
Theta = 0
Camera = np.zeros(3)

def MultiplyMatrixVector(vec3di, mat4x4):
    vec3d = vec3di.copy()
    vec3d[0] = vec3di[0] * mat4x4[0][0] + vec3di[1] * mat4x4[1][0] + vec3di[2] * mat4x4[2][0] + mat4x4[3][0]
    vec3d[1] = vec3di[0] * mat4x4[0][1] + vec3di[1] * mat4x4[1][1] + vec3di[2] * mat4x4[2][1] + mat4x4[3][1]
    vec3d[2] = vec3di[0] * mat4x4[0][2] + vec3di[1] * mat4x4[1][2] + vec3di[2] * mat4x4[2][2] + mat4x4[3][2]
    w = vec3di[0] * mat4x4[0][3] + vec3di[1] * mat4x4[1][3] + vec3di[2] * mat4x4[2][3] + mat4x4[3][3]
    if w != 0:
        vec3d[0] /= w;        vec3d[1] /= w;        vec3d[2] /= w
    
    return vec3d

def read_obj(fileName):
    vertices = [[0,0,0]]
    mesh = []
    f = open(fileName)
    for line in f:
        if line[:2] == "v ":
            index1 = line.find(" ") + 1
            index2 = line.find(" ", index1 + 1)
            index3 = line.find(" ", index2 + 1)

            vertex = [float(line[index1:index2]), float(line[index2:index3]), float(line[index3:-1])]
            vertices.append(vertex)

        elif line[0] == "f":
            index1 = line.find(" ") + 1
            index2 = line.find(" ", index1 + 1)
            index3 = line.find(" ", index2 + 1)

            face = (int(line[index1:index2]), int(line[index2:index3]), int(line[index3:-1]))
            mesh.append([vertices[face[0]], vertices[face[1]], vertices[face[2]]])

    f.close()

    return mesh


# meshCube = [

# #  SOUTH
# [ [0, 0, 0],    [0, 1, 0],    [1, 1, 0 ]],
# [ [0, 0, 0],    [1, 1, 0],    [1, 0, 0 ]],

# #  EAST                                                      
# [ [1, 0, 0],    [1, 1, 0],    [1, 1, 1 ]],
# [ [1, 0, 0],    [1, 1, 1],    [1, 0, 1 ]],

# #  NORTH                                                     
# [ [1, 0, 1],    [1, 1, 1],    [0, 1, 1 ]],
# [ [1, 0, 1],    [0, 1, 1],    [0, 0, 1 ]],

# #  WEST                                                      
# [ [0, 0, 1],    [0, 1, 1],    [0, 1, 0 ]],
# [ [0, 0, 1],    [0, 1, 0],    [0, 0, 0 ]],

# #  TOP                                                       
# [ [0, 1, 0],    [0, 1, 1],    [1, 1, 1 ]],
# [ [0, 1, 0],    [1, 1, 1],    [1, 1, 0 ]],

# #  BOTTOM                                                    
# [ [1, 0, 1],    [0, 0, 1],    [0, 0, 0 ]],
# [ [1, 0, 1],    [0, 0, 0],    [1, 0, 0 ]],

# ]
meshCube = read_obj('VideoShip.obj')

#  Projection Matrix
Near = 0.1
Far = 1000
Fov = 90
AspectRatio = SCREEN_HEIGHT / SCREEN_WIDTH
FovRad = 1 / np.tan(Fov * 0.5 / 180 * 3.14159)
matProj = np.zeros((4,4))
matProj[0][0] = AspectRatio * FovRad
matProj[1][1] = FovRad
matProj[2][2] = Far / (Far - Near)
matProj[3][2] = (-Far * Near) / (Far - Near)
matProj[2][3] = 1
matProj[3][3] = 0

while running:
    elapsed_time = clock.tick()/1000
    for event in pg.event.get():
        if event.type == pg.QUIT:
            running = False
        if event.type == pg.KEYDOWN:
            if event.key == pg.K_ESCAPE:
                running = False
    
    #  Clear Screen
    surf.fill([0,0,0])

    #  Set up rotation matrices
    Theta += 1 * elapsed_time

    #  Rotation Z
    matRotZ = np.zeros((4,4))
    matRotZ[0][0] = np.cos(Theta)
    matRotZ[0][1] = np.sin(Theta)
    matRotZ[1][0] = -np.sin(Theta)
    matRotZ[1][1] = np.cos(Theta)
    matRotZ[2][2] = 1
    matRotZ[3][3] = 1

    #  Rotation X
    matRotX = np.zeros((4,4))
    matRotX[0][0] = 1
    matRotX[1][1] = np.cos(Theta * 0.5)
    matRotX[1][2] = np.sin(Theta * 0.5)
    matRotX[2][1] = -np.sin(Theta * 0.5)
    matRotX[2][2] = np.cos(Theta * 0.5)
    matRotX[3][3] = 1
    
    #  Draw Triangles
    TrianglesToRaster = []
    z_order = []
    for triangle in meshCube:
        triProjected = triangle.copy()
        triTranslated = triangle.copy()
        triRotatedZ = triangle.copy()
        triRotatedZX = triangle.copy()

        #  Rotate in Z-Axis
        triRotatedZ[0] = MultiplyMatrixVector(triangle[0], matRotZ)
        triRotatedZ[1] = MultiplyMatrixVector(triangle[1], matRotZ)
        triRotatedZ[2] = MultiplyMatrixVector(triangle[2], matRotZ)

        #  Rotate in X-Axis
        triRotatedZX[0] = MultiplyMatrixVector(triRotatedZ[0], matRotX)
        triRotatedZX[1] = MultiplyMatrixVector(triRotatedZ[1], matRotX)
        triRotatedZX[2] = MultiplyMatrixVector(triRotatedZ[2], matRotX)


        #  Offset into the screen
        triTranslated = triRotatedZX.copy()
        triTranslated[0][2] = triRotatedZX[0][2] + 8
        triTranslated[1][2] = triRotatedZX[1][2] + 8
        triTranslated[2][2] = triRotatedZX[2][2] + 8
        # Use Cross-Product to get surface normal

        normal, line1, line2 = np.zeros(3), np.zeros(3), np.zeros(3)
        line1[0] = triTranslated[1][0] - triTranslated[0][0]
        line1[1] = triTranslated[1][1] - triTranslated[0][1]
        line1[2] = triTranslated[1][2] - triTranslated[0][2]

        line2[0] = triTranslated[2][0] - triTranslated[0][0]
        line2[1] = triTranslated[2][1] - triTranslated[0][1]
        line2[2] = triTranslated[2][2] - triTranslated[0][2]

        normal[0] = line1[1] * line2[2] - line1[2] * line2[1]
        normal[1] = line1[2] * line2[0] - line1[0] * line2[2]
        normal[2] = line1[0] * line2[1] - line1[1] * line2[0]

        # It's normally normal to normalise the normal
        lenght = np.sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]) +1e-16
        normal[0] /= lenght; normal[1] /= lenght; normal[2] /= lenght

        # if (normal.z < 0)
        if (normal[0] * (triTranslated[0][0] - Camera[0]) + 
            normal[1] * (triTranslated[0][1] - Camera[1]) +
            normal[2] * (triTranslated[0][2] - Camera[2]) < 0.0):
            
            # Illumination
            light_direction = [0, 0, -1]
            lenght = np.sqrt(light_direction[0]**2 + light_direction[1]**2 + light_direction[2]**2)
            light_direction[0] /= lenght; light_direction[1] /= lenght; light_direction[2] /= lenght

            # How similar is normal to light direction
            dp = normal[0] * light_direction[0] + normal[1] * light_direction[1] + normal[2] * light_direction[2]
            color = np.abs(np.ones(3)*dp*255)
            #  Project triangles from 3D --> 2D
            triProjected[0] = MultiplyMatrixVector(triTranslated[0], matProj)
            triProjected[1] = MultiplyMatrixVector(triTranslated[1], matProj)
            triProjected[2] = MultiplyMatrixVector(triTranslated[2], matProj)

            #  Scale into view
            triProjected[0][0] += 1; triProjected[0][1] += 1
            triProjected[1][0] += 1; triProjected[1][1] += 1
            triProjected[2][0] += 1; triProjected[2][1] += 1
            triProjected[0][0] *= 0.5 * SCREEN_WIDTH
            triProjected[0][1] *= 0.5 * SCREEN_HEIGHT
            triProjected[1][0] *= 0.5 * SCREEN_WIDTH
            triProjected[1][1] *= 0.5 * SCREEN_HEIGHT
            triProjected[2][0] *= 0.5 * SCREEN_WIDTH
            triProjected[2][1] *= 0.5 * SCREEN_HEIGHT

            #  Rasterize triangle
            point0 = [triProjected[0][0], triProjected[0][1]]
            point1 = [triProjected[1][0], triProjected[1][1]]
            point2 = [triProjected[2][0], triProjected[2][1]]
            z_order.append(-(triProjected[0][2] + triProjected[1][2] + triProjected[2][2])/3)
            TrianglesToRaster.append([point0, point1, point2, color])
            #pg.draw.polygon(surf, color, [point0, point1, point2])
            #pg.draw.lines(surf, [0, 0, 255], 1, [point0, point1, point2])
    
    for index in np.argsort(np.asarray(z_order),):
        triangle = TrianglesToRaster[index]
        pg.draw.polygon(surf, triangle[3], [triangle[0], triangle[1], triangle[2]])
    screen.blit(surf, (0,0))
    pg.display.update()
    pg.display.set_caption(str(round(1/(elapsed_time+1e-16), 1)))

pg.quit()
