import numpy as np

def closest_point_on_line(a, b, p):
    """ Find the closest point on the line segment defined by points A and B to point P """
    ap = p - a
    ab = b - a
    result = a + np.dot(ap, ab) / np.dot(ab, ab) * ab
    return result

def collision_detector(line_point1, line_point2, sphere_center, sphere_radius):
    """ Check if the line defined by two points collides with the sphere """
    closest_point = closest_point_on_line(np.array(line_point1), np.array(line_point2), np.array(sphere_center))
    distance = np.linalg.norm(closest_point - np.array(sphere_center))
    
    #check if the sphere is beyond the line segment
    if np.dot(closest_point - np.array(line_point1), closest_point - np.array(line_point2)) > 0:
        return False
    else:
        
         return distance <= sphere_radius

# # Example usage:
# line_point1 = [-3, 0, 0]
# line_point2 = [-0.5, 0, 0]
# sphere_center = [0, 1, 0]
# sphere_radius = 1

# collision = collision_detector(line_point1, line_point2, sphere_center, sphere_radius)
# print("Collision:", collision)
