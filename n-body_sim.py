import matplotlib.pyplot as plt
import numpy as np
import copy

# Define constants
GRAV_CONSTANT = 6.67408e-11

class Body(object):
    """ Defines a body in 3D space with mass m, initial velocity vector v and radius r. All units are metric. """
    def __init__(self, name, position, mass, velocity, radius, color='b'):
        self.name = name
        self.pos = position.astype(np.float64)
        self.m = mass
        self.v = velocity.astype(np.float64)
        self.r = float(radius)
        self.color = color


    def __repr__(self):
        return 'Name: ' + self.name + ' position: ' + str(self.pos) + ' velocity:' + str(self.v)

    def __str__(self):
        """ What will be printed when this object is printed. """

        return 'Name: ' + self.name + ' position: ' + str(self.pos) + ' velocity:' + str(self.v)



def getForce(obj1, obj2, eps=1e-12):
    """ Calculates the force between two bodies. """

    # Calculate the distance between two positions ("eps" saves the day when 2 objects are on the exactly same 
    # spot)
    dist = np.linalg.norm(obj1.pos - obj2.pos) + eps

    # Calculate the force
    force_vect = -GRAV_CONSTANT*obj1.m*obj2.m * (obj1.pos - obj2.pos)/dist**3

    return force_vect


def timeStep(obj_list, time_step, collisions=False):
    """ Make one time step with the given objects. The time step is given in seconds. """

    # Iterate over all bodies to update the position
    for i in range(len(obj_list)):
        obj_list[i].pos += obj_list[i].v*time_step/2

    # Iterate over all bodies and calculate new velocity vectors
    for i in range(len(obj_list)):

        obj1 = obj_list[i]

        force_sum = np.zeros(3)
        
        for j in range(len(obj_list)):

            obj2 = obj_list[j]

            # Do not calculte for the same object
            if i != j:
                force_sum += getForce(obj1, obj2)

        # Convert force to velocity
        a = force_sum/obj1.m
        v_diff = a * time_step

        # Update velocity
        obj_list[i].v += v_diff

    # Update object position
    for i in range(len(obj_list)):
        obj_list[i].pos += obj_list[i].v*time_step/2


    if collisions:
        # Check for colisions
        for i, obj1 in enumerate(obj_list):

            for j, obj2 in enumerate(obj_list):
                if obj1.name != obj2.name:

                    # Check the distance between bodies
                    dist = np.linalg.norm(obj1.pos - obj2.pos)

                    # Merge bodies if they collide
                    if (dist <= obj1.r or dist <= obj2.r):

                        # Sum the masses
                        mass = obj1.m + obj2.m

                        # Sum the momenta
                        momentum = obj1.m*obj1.v + obj2.m*obj2.v

                        # Calculate the final velocity
                        velocity = momentum/mass

                        # Calculate the new position
                        position = (obj1.m*obj1.pos + obj2.m*obj2.pos)/mass

                        # Calculate new object radius (assume they are spheres)
                        r = np.power(obj1.r**3 + obj2.r**3, 1.0/3)

                        # Create the new merged body
                        new_obj = Body(obj1.name+obj2.name, position, mass, velocity, r)


                        # Remove old bodies
                        obj_list.pop(i)

                        if i <= j:
                            j -= 1
                        obj_list.pop(j)

                        # Add new body to the list
                        obj_list.append(new_obj)

                        break


    return obj_list


def plotBodies(obj_list):

    # Clear the plot
    plt.clf()

    # Set plot limits
    
    # plt.xlim((-2.5e11, 2.5e11))
    # plt.ylim((-2.5e11, 2.5e11))

    plt.xlim((-20, 20))
    plt.ylim((-20, 20))

    # Plot object positions
    for obj in obj_list:

        # Plot position
        plt.scatter(obj.pos[0], obj.pos[1], s=norm_param*np.pi*obj.r**2, c=obj.color)

        v_scale = 20

        # Plot velocity vector
        plt.arrow(obj.pos[0], obj.pos[1], v_scale*obj.v[0], v_scale*obj.v[1], head_width=0.03*v_scale, head_length=0.05*v_scale, color='r')

    # Update the plot
    plt.draw()
    plt.pause(0.0001)



if __name__ == '__main__':

    # Testing

    obj1 = Body('1', np.array([0, 0, 0]), 1e8, np.array([0, 0.12, 0]), 2)
    obj2 = Body('2', np.array([10, 10, 0]), 1e7, np.array([0, 0, 0]), 1)
    obj3 = Body('3', np.array([5, 0, 0]), 1e9, np.array([0.01, 0, 0]), 3)
    obj4 = Body('4', np.array([-10, -10, 0]), 1e8, np.array([0.1, 0.1, 0]), 2)
    obj5 = Body('5', np.array([0, -10, 0]), 1e7, np.array([0.1, 0, 0]), 1)

    obj_list = [obj1, obj2, obj3, obj4, obj5]


    # # The Solar System
    # sun     = Body('sun',     np.array([0,       0,       0]), 2e30,   np.array([0,        0,       0]), 10, color='y')
    # mercury = Body('mercury', np.array([0,       5.8e10,  0]), 3.3e23, np.array([-47.4e3,  0,       0]), 2,  color='orange')
    # venus   = Body('venus',   np.array([108.2e9, 0,       0]), 4.9e24, np.array([0,        35e3,    0]), 3,  color='y')
    # earth   = Body('earth',   np.array([150e9,   0,       0]), 6e24,   np.array([0,        29865.3, 0]), 5,  color='b')
    # mars    = Body('mars',    np.array([0,       227.9e9, 0]), 6.4e23, np.array([-24130.8, 0,       0]), 4,  color='r')

    # obj_list = [sun, mercury, venus, earth, mars]

    # Find normalizing parameter for the object size plotting
    max_radius = max([obj.r for obj in obj_list])
    norm_param = (20.0 / max_radius)**2

    # Prepare for plotting in real time
    plt.ion()
    plt.show()

    plotBodies(obj_list)

    for i in range(500):
        obj_list = timeStep(obj_list, 0.5, collisions=True)

        # Change the number after % to reduce the plotting 'framerate'
        if i%5 == 0:
            plotBodies(obj_list)