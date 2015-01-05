module MixedVoxels

export intersection_points, polygon_area

# intersect a plane with the einheits voxel
function intersection_points(n_p, d)
	
	intersectionpoints::Vector(Float64(3)) = []

	# minimal distance of an intersection point from a cube corner
	eps = 1e-7

	# no intersection possible if
	if abs(d) > sqrt(3.0)
		return intersectionpoints
	end

	# consider each corner (x,y,z) in {0,1}^3
	for z in (0, 1)
		for y in (0, 1)
			for x in (0, 1)
				# corner
				c = [x, y, z]

				# check if corner is in plane and plane inside the voxel
				if # TODO: wie equalszero()?
				end

				# skip corner (1,1,1)
				if c == [1, 1, 1] continue end

				# 3 possible edges for this corner
				for j = 1:3

					# edge direction
					e = [0, 0, 0]

					# only consider edges in positive directions
					if c[j] == 1
						continue
					else
						e[j] = 1
					end

					# scalar product of plane normal and edge direction
					denom = dot(n_p, e)

					# check edge for intersection with plane
					if denom != 0
						# TODO: wie equalszero()?
						lambda = (d - dot(n_p, c)) / denom
					else
						# edge coplanar with plane
						lambda = -1.0
					end

					# test for intersection within cube
					if (lambda > eps) && (lambda < 1.0-eps)
						
						# intersection
						point = c + lambda * e

						# add intersection point
						push!(intersectionpoints, point)

						# there's six points of intersection at maximum
						if length(intersectionpoints) == 6
							return intersectionpoints
						end
					end

				end # continue with another edge
			end
		end
	end # continue with another corner

	return intersectionpoints
end

# sort vertices of a polygon in positive (counterclockwise) direction
# around the normal vector n_p
function sort()

end

# Gauss's shoelace formula for the area of a polygon
# points must be ordered in positive direction (counterclockwise)
# use sort() for that
function polygon_area(vertices::Vector(Float64(3)))

	for i in 1:length(vertices)
		v += cross(vertices[i], vertices[mod(i+1, length(vertices))])
	end

	return 0.5 * norm(v)
end

end