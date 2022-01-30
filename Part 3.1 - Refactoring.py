import pygame as pg
import numpy as np

def Matrix_MakeRotationX(AngleRad):
    matrix = np.identity(4)
    matrix[1][1] = np.cos(AngleRad)
    matrix[1][2] = np.sin(AngleRad)
    matrix[2][1] = -np.sin(AngleRad)
    matrix[2][2] = np.cos(AngleRad)

    return matrix
	
def Matrix_MakeRotationY(AngleRad):
    matrix = np.identity(4)
    matrix[0][0] = np.cos(AngleRad)
    matrix[0][2] = np.sin(AngleRad)
    matrix[2][0] = -np.sin(AngleRad)
    matrix[2][2] = np.cos(AngleRad)
	
    return matrix

def Matrix_MakeRotationZ(AngleRad):
    matrix = np.identity(4)
    matrix[0][0] = np.cos(AngleRad)
    matrix[0][1] = np.sin(AngleRad)
    matrix[1][0] = -np.sin(AngleRad)
    matrix[1][1] = np.cos(AngleRad)

    return matrix

def Matrix_MakeTranslation(x, y, z):
    matrix = np.identity(4)
    matrix[3][0] = x
    matrix[3][1] = y
    matrix[3][2] = z
	
    return matrix

def Matrix_MakeProjection(FovDegrees, AspectRatio, Near, Far):
    FovRad = 1 / np.tan(FovDegrees * 0.5 / 180 * 3.14159)
    matrix = np.identity(4)
    matrix[0][0] = AspectRatio * FovRad
    matrix[1][1] = FovRad
    matrix[2][2] = Far / (Far - Near)
    matrix[3][2] = (-Far * Near) / (Far - Near)
    matrix[2][3] = 1
    matrix[3][3] = 0
	
    return matrix

def Vector_Normalise(vec):
    lenght = np.sqrt(vec[0]**2 + vec[1]**2 + vec[1]**2)
    if lenght > 0:
        vec = vec/lenght

    return vec

def Matrix_PointAt(pos, target, up):
    # Calculate new forward direction
    newForward = target - pos
    newForward = Vector_Normalise(newForward)

    # Calculate new Up direction
    a = np.multiply(newForward, np.dot(up, newForward))
    newUp = up - a
    newUp = Vector_Normalise(newUp)
    # New Right direction is easy, its just cross product
    newRight = np.cross(newUp[:3], newForward[:3])

    # Construct Dimensioning and Translation Matrix	
    matrix = np.zeros((4,4))
    matrix[0][0] = newRight[0];	    matrix[0][1] = newRight[1];     matrix[0][2] = newRight[2];	matrix[0][3] = 0
    matrix[1][0] = newUp[0];		matrix[1][1] = newUp[1];		matrix[1][2] = newUp[2];		matrix[1][3] = 0
    matrix[2][0] = newForward[0];	matrix[2][1] = newForward[1];	matrix[2][2] = newForward[2];	matrix[2][3] = 0
    matrix[3][0] = pos[0];		    matrix[3][1] = pos[1];		    matrix[3][2] = pos[2];		matrix[3][3] = 1
    
    return matrix


def read_obj(fileName):
    vertices = []
    mesh = []
    f = open(fileName)
    for line in f:
        if line[:2] == "v ":
            index1 = line.find(" ") + 1
            index2 = line.find(" ", index1 + 1)
            index3 = line.find(" ", index2 + 1)

            vertex = [float(line[index1:index2]), float(line[index2:index3]), float(line[index3:-1]), 1]
            vertices.append(np.asarray(vertex))

        elif line[0] == "f":
            index1 = line.find(" ") + 1
            index2 = line.find(" ", index1 + 1)
            index3 = line.find(" ", index2 + 1)

            face = [int(line[index1:index2]) - 1, int(line[index2:index3]) - 1, int(line[index3:-1]) - 1]
            mesh.append([vertices[face[0]], vertices[face[1]], vertices[face[2]]])

    f.close()

    return mesh


SCREEN_WIDTH, SCREEN_HEIGHT = 800, 600

pg.init()
screen = pg.display.set_mode((SCREEN_WIDTH,SCREEN_HEIGHT))
running = True
clock = pg.time.Clock()
surf = pg.surfarray.make_surface(np.zeros((SCREEN_WIDTH,SCREEN_HEIGHT, 3)))
Theta = 0
Camera = np.zeros(4)
Near = 0.1
Far = 1000
FovDegrees = 90
AspectRatio = SCREEN_HEIGHT / SCREEN_WIDTH
meshCube = read_obj('VideoShip.obj')
matProj = Matrix_MakeProjection(FovDegrees, AspectRatio, Near, Far)

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
    Theta += 1 * elapsed_time # Uncomment to spin me right round baby right round
	
    # Set up "World Tranmsform"
    matRotZ = Matrix_MakeRotationZ(Theta * 0.5)
    matRotX = Matrix_MakeRotationX(Theta)

    matTrans = Matrix_MakeTranslation(0, 0, 10)

    # Form World Matrix
    matWorld = np.matmul(matRotZ, matRotX) # Transform by rotation
    matWorld = np.matmul(matWorld, matTrans) # Transform by translation

    #  Draw Triangles
    TrianglesToRaster = []
    z_order = []
    for triangle in meshCube:
        triTransformed = triangle.copy()

        #  Rotate in Z-Axis
        triTransformed[0] = np.matmul(triangle[0], matWorld)
        triTransformed[1] = np.matmul(triangle[1], matWorld)
        triTransformed[2] = np.matmul(triangle[2], matWorld)

        
        # Use Cross-Product to get surface normal
        line1 = triTransformed[1] - triTransformed[0]
        line2 = triTransformed[2] - triTransformed[0]

        # It's normally normal to normalise the normal
        normal = Vector_Normalise(np.cross(line1[:3], line2[:3]))

        CameraRay = triTransformed[0] - Camera

        # if (normal.z < 0)
        if (np.dot(normal, CameraRay[:3]) < 0):
            
            # Illumination
            light_direction = np.asarray([0, 0, -1])
            light_direction = Vector_Normalise(light_direction)
            
            # How similar is normal to light direction
            dp = min(1, max(0.1, np.dot(light_direction, normal)))
            color = np.ones(3)*dp*255

            # Project triangles from 3D --> 2D
            triProjected = triTransformed.copy()
            triProjected[0] = np.matmul(triTransformed[0], matProj)
            triProjected[1] = np.matmul(triTransformed[1], matProj)
            triProjected[2] = np.matmul(triTransformed[2], matProj)

            # Scale into view, we moved the normalising into cartesian space
			      # out of the matrix.vector function from the previous videos, so
			      # do this manually
            triProjected[0] = triProjected[0]/triProjected[0][3]
            triProjected[1] = triProjected[1]/triProjected[1][3]
            triProjected[2] = triProjected[2]/triProjected[2][3]

            # X/Y are inverted so put them back
            triProjected[0][0] *= -1
            triProjected[1][0] *= -1
            triProjected[2][0] *= -1
            triProjected[0][1] *= -1
            triProjected[1][1] *= -1
            triProjected[2][1] *= -1
            
            #  Scale into view
            OffsetView = np.asarray([1,1,0, 0])
            triProjected[0] += OffsetView
            triProjected[1] += OffsetView
            triProjected[2] += OffsetView
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
    
    for index in np.argsort(np.asarray(z_order)):
        triangle = TrianglesToRaster[index]
        pg.draw.polygon(surf, triangle[3], [triangle[0], triangle[1], triangle[2]])

    screen.blit(surf, (0,0))
    pg.display.update()
    pg.display.set_caption(str(round(1/(elapsed_time+1e-16), 1)))

pg.quit()
