module Archimedes

using Unitful
using Luxor

export Point2D,
       PointMass,
       mass,
       getx,
       gety,
       radius,
       bbox,
       centroid,
       gravitycenter,
       draw,
       drawing,
       show_drawing,
       BoxShip,
       corners,
       labeled_point,
       buoyancycenter,
       toship

immutable Point2D
    x::typeof(1.0u"m")
    y::typeof(1.0u"m")
end
getx(p::Point2D) = p.x
gety(p::Point2D) = p.y

immutable PointMass
    coordinates::Point2D
    mass::typeof(1.0u"kg")
    radius::typeof(1.0u"m")
end
PointMass(x,y,m,r) = PointMass(Point2D(x,y),m,r)
coords(p::PointMass) = p.coordinates
mass(p::PointMass) = p.mass
getx(p::PointMass) = getx(coords(p))
gety(p::PointMass) = gety(coords(p))
radius(p::PointMass) = p.radius

"""
Get the bounding box (xmin, xmax, ymin, ymax) of an array of pointmasses
"""
bbox(pointmasses) = (extrema(getx.(pointmasses))..., extrema(gety.(pointmasses))...)

"""
Get the centroid of an array of pointmasses
"""
centroid(pointmasses) = Point2D(mean(getx.(pointmasses)), mean(gety.(pointmasses)))

"""
Get the center of gravity
"""
function gravitycenter(pointmasses)
  m = sum(mass.(pointmasses))
  x = sum(getx.(pointmasses) .* mass.(pointmasses)) / m
  y = sum(gety.(pointmasses) .* mass.(pointmasses)) / m
  r = sqrt(sum(radius.(pointmasses) .^ 2))
  return PointMass(x, y, m, r)
end

"""
Helper to map the points with units to pixel-coordinates
"""
immutable CoordMapping
  origin::Point2D
  scaling::typeof(1.0/(1.0u"m"))
end

function CoordMapping(width::Number, bbox::Tuple)
  (xmin,xmax,ymin,ymax) = bbox
  center = Point2D(mean((xmin,xmax)), mean((ymin,ymax)))
  scaling = width/(xmax-xmin)
  return CoordMapping(center, scaling)
end

scale(x, mapping::CoordMapping) = x*mapping.scaling

function remap(p, mapping)
  (xoff, yoff) = (getx(mapping.origin), gety(mapping.origin))
  return Luxor.Point(((getx(p)-xoff, -(gety(p)-yoff)).*mapping.scaling)...)
end

# filename for the SVG temp file
const tempfile = "_archimedes_drawing.svg"

"""
Initialize a drawing, returing the point mapping
"""
function drawing(bbox, figwidth=300)
  (xmin,xmax,ymin,ymax) = bbox
  ar = (xmax-xmin)/(ymax-ymin)
  figheight = figwidth/ar
  margin = 30
  bottom_padding = 40
  d = Drawing(figwidth+2*margin, figheight+2*margin+bottom_padding, tempfile)
  origin()
  background("white")
  return CoordMapping(figwidth, bbox)
end

"""
Draws a pointmass into the current context
"""
function draw(p::PointMass, mapping::CoordMapping; hue="black")
  sethue(hue)
  c = remap(p,mapping)
  r = scale(radius(p), mapping)
  circle(c, r, :fill)
  spacing = 5.0
  step = 10.0
  text("x = $(getx(p))", Point(c.x,c.y+r+spacing), halign=:center, valign=:top)
  text("y = $(gety(p))", Point(c.x,c.y+r+spacing+step), halign=:center, valign=:top)
  text("m = $(mass(p))", Point(c.x,c.y+r+spacing+2*step), halign=:center, valign=:top)
end

"""
Shows the drawing and removes the temp file
"""
function show_drawing()
  finish()
  preview()
  rm(tempfile)
end

immutable BoxShip
  width::typeof(1.0u"m")
  height::typeof(1.0u"m")
  draft::typeof(1.0u"m")
  KG::typeof(1.0u"m")
  heel::Float64
end

rotmat(s::BoxShip) = [cos(s.heel) -sin(s.heel); sin(s.heel) cos(s.heel)]

"""
Convert position p to ship coordinates
"""
function toship(p, s)
  rmat = rotmat(s)
  rmat[1,2] = -rmat[1,2]
  rmat[2,1] = -rmat[2,1]
  return Point2D((rmat*[getx(p), gety(p)])...)
end

function corners(s::BoxShip)
  rotm = rotmat(s)
  x = s.width/2
  ymax = s.height-s.draft
  ymin = -s.draft
  points = rotm*[-x x x -x; ymin ymin ymax ymax]
  return Point2D.(points[1,:], points[2,:])
end

function waterline(s::BoxShip)
  x = 1.4*s.width/2
  y = 0.0u"m"
  return (Point2D(-x,y), Point2D(x,y))
end

function gravitycenter(s::BoxShip)
  g_coords = rotmat(s)*[0.0u"m", s.KG-s.draft]
  return Point2D(g_coords...)
end

bbox(s::BoxShip) = 1.4 .* (-s.width/2, s.width/2, -s.draft, s.height-s.draft)

"""
Corners of the underwater part
"""
function carene(s::BoxShip)
  ship_pts = corners(s)
  x1 = getx(ship_pts[4])
  y1 = gety(ship_pts[4])
  x2 = getx(ship_pts[1])
  y2 = gety(ship_pts[1])
  r = (x2-x1)/(y2-y1)
  wl_left = Point2D(-y1*r+x1, 0.0u"m")
  wl_right = Point2D(-gety(ship_pts[3])*r+getx(ship_pts[3]), 0.0u"m")
  return (ship_pts[1], ship_pts[2], wl_right, wl_left)
end

"""
Make a 3-vector out of a Point2D
"""
vec3(p::Point2D) = [getx(p), gety(p), 0.0u"m"]

"""
Area of a triangle
"""
function area(pts)
  v = vec3.(pts)
  return abs(cross(v[3]-v[1], v[2]-v[1])[3])/2
end

function buoyancycenter(s::BoxShip)
  trap = carene(s)
  c = centroid(trap)
  triags = []
  N = length(trap)
  for i in linearindices(trap)
    push!(triags, [c, trap[i], trap[i%N+1]])
  end
  areas = area.(triags)
  centroids = centroid.(triags)
  A = sum(areas)
  x = sum(getx.(centroids) .* areas) / A
  y = sum(gety.(centroids) .* areas) / A
  return Point2D(x,y)
end


function labeled_point(p, label, hue="black", radius=5.0)
  sethue(hue)
  circle(p, radius, :fill)
  text(label, Point(p.x+1.2*radius,p.y+1.2*radius), halign=:center, valign=:top)
end

function draw(s::BoxShip, figwidth=300)
  ship_pts = corners(s)
  m = drawing(bbox(s), figwidth)
  p = remap.(ship_pts, m)
  sethue("black")
  line(p[1], p[2], :stroke)
  line(p[2], p[3], :stroke)
  line(p[4], p[1], :stroke)
  labeled_point(remap(gravitycenter(s),m), "G")
  labeled_point(remap(buoyancycenter(s),m), "B", "red")
  wl = remap.(waterline(s), m)
  sethue("blue")
  line(wl[1], wl[2], :stroke)
  return m
end

end # module
