

from scipy.spatial.transform import Rotation as R
import numpy as np
    
def translate_and_rotate_to_x_axis(P1, P2):
    """
    Translate the line defined by points P1 and P2 to the origin and rotate it to align with the X-axis.
    """
    # Translate the points
    P1_translated = np.array([0, 0, 0])
    P2_translated = np.array(P2) - np.array(P1)

    # Direction vector of the line
    direction_vector = P2_translated / np.linalg.norm(P2_translated)

    # Calculate the angles with the X, Y, Z axes
    x_axis = np.array([1, 0, 0])
    angle_with_x_axis = np.arccos(np.dot(direction_vector, x_axis))

    # Rotation axis (cross product with x-axis)
    rotation_axis = np.cross(direction_vector, x_axis)
    rotation_axis = rotation_axis / np.linalg.norm(rotation_axis)

    # Create rotation object
    rotation = R.from_rotvec(rotation_axis * angle_with_x_axis)

    # Rotate the line
    P1_rotated = rotation.apply(P1_translated)
    P2_rotated = rotation.apply(P2_translated)

    return P1_rotated, P2_rotated

def translate_and_rotate_point_to_x_axis(P, P1, P2):
    """
    Apply the same transformation used in translate_and_rotate_to_x_axis to a single point P.
    """
    # Translate the point
    P1_translated = np.array([0, 0, 0])
    P_translated = np.array(P) - np.array(P1)

    # Reuse the calculation of direction vector, angle, and rotation from the line function
    # Direction vector of the line
    P2_translated = np.array(P2) - np.array(P1)
    direction_vector = P2_translated / np.linalg.norm(P2_translated)

    # Calculate the angles with the X, Y, Z axes
    x_axis = np.array([1, 0, 0])

    # Check if direction vector is parallel to the x-axis
    if np.allclose(direction_vector, x_axis) or np.allclose(direction_vector, -x_axis):
        # No rotation needed as the line is already aligned with the x-axis
        return P_translated

    angle_with_x_axis = np.arccos(np.dot(direction_vector, x_axis))

    # Rotation axis (cross product with x-axis)
    rotation_axis = np.cross(direction_vector, x_axis)

    # Check if rotation axis is zero vector
    if not np.any(rotation_axis):
        # This means the line is parallel to the x-axis, and no rotation is needed
        return P_translated

    rotation_axis = rotation_axis / np.linalg.norm(rotation_axis)

    # Create rotation object
    rotation = R.from_rotvec(rotation_axis * angle_with_x_axis)
    P1_rotated = rotation.apply(P1_translated)
    P2_rotated = rotation.apply(P2_translated)

    # Rotate the point
    P_rotated = rotation.apply(P_translated)
    P1_rotated = rotation.apply(P1_translated)
    P2_rotated = rotation.apply(P2_translated)
    return P_rotated, P1_rotated, P2_rotated

def collision_detector(P1,P2, P, P_radius):
    P_rotated, P1_rotated, P2_rotated = translate_and_rotate_point_to_x_axis(P1, P2, P)
    # print(P_rotated)
    # print(P1_rotated)
    # print(P2_rotated)
    # print("Distance to line: ", np.linalg.norm(P_rotated[1:]))
    if (P_rotated[0] < 0) or (P_rotated[0] > P2_rotated[0]): 
        return False
    if np.linalg.norm(P_rotated[1:]) <= P_radius: #here we check the y and z coordinates of the 
        return True
    else:
        return False
    

#uncomment the following lines to test the functions

# # Example usage
# P = [0, 4,1]  # Example point
# P1 = [0, 1, 0]  # Line start
# P2 = [0, 3, 0]  # Line end
# transformed_point = translate_and_rotate_point_to_x_axis(P, P1, P2)
# transformed_P1, transformed_P2 = translate_and_rotate_to_x_axis(P1, P2)

# print(collision_detector(P, P2, P1, 1))

# #plotting
# import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D
# import numpy as np

# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')

# # Plot the line
# ax.plot([P1[0], P2[0]], [P1[1], P2[1]], [P1[2], P2[2]], color='k')

# # Plot the point
# ax.scatter(P[0], P[1], P[2], color='r')

# # Plot the transformed line

# ax.plot([transformed_P1[0], transformed_P2[0]], [transformed_P1[1], transformed_P2[1]], [transformed_P1[2], transformed_P2[2]], color='b')

# # Plot the transformed point

# ax.scatter(transformed_point[0], transformed_point[1], transformed_point[2], color='g')

# plt.savefig("transformation.png")
