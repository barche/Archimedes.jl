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
       toship,
       isocarene,
       metacenter,
       setopacity,
       puttext

immutable Point2D
    x::typeof(1.0u"m")
    y::typeof(1.0u"m")
end
getx(p::Point2D) = p.x
gety(p::Point2D) = p.y
Base.:-(a::Point2D, b::Point2D) = Point2D(a.x-b.x, a.y-b.y)
Base.:+(a::Point2D, b::Point2D) = Point2D(a.x+b.x, a.y+b.y)
Base.:*(a::Point2D, b::Number) = Point2D(a.x*b, a.y*b)
Base.:*(b::Number, a::Point2D) = Point2D(a.x*b, a.y*b)
Base.:/(a::Point2D, b::Number) = Point2D(a.x/b, a.y/b)
Base.norm(a::Point2D) = sqrt(a.x^2 + a.y^2)

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

"""
Initialize a drawing, returing the point mapping
"""
function drawing(bbox, figwidth=300)
  (xmin,xmax,ymin,ymax) = bbox
  ar = (xmax-xmin)/(ymax-ymin)
  figheight = figwidth/ar
  margin = 30
  bottom_padding = 40
  d = Drawing(figwidth+2*margin, figheight+2*margin+bottom_padding, :svg)
  origin()
  background("white")
  return CoordMapping(figwidth, bbox)
end

function puttext(txt, p::Point2D, mapping::CoordMapping)
  sethue("black")
  c = remap(p,mapping)
  text(txt, c, halign=:center, valign=:top)
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

function draw(p::Point2D, mapping::CoordMapping, r=5.0; hue="black", label="")
  sethue(hue)
  c = remap(p,mapping)
  circle(c, r, :fill)
  if label != ""
    text(label, Point(c.x+1.2*r,c.y+1.2*r), halign=:center, valign=:top)
  end
end

"""
Shows the drawing and removes the temp file
"""
function show_drawing()
  finish()
  return Luxor.currentdrawing
end

immutable BoxShip
  width::typeof(1.0u"m")
  height::typeof(1.0u"m")
  draft::typeof(1.0u"m")
  KG::typeof(1.0u"m")
  heel::Float64
  vshift::typeof(1.0u"m")

  BoxShip(width, height, draft, KG, heel=0.0, vshift=0.0u"m") = new(width,height,draft,KG,heel,vshift)
end

forward_trans(p, s::BoxShip) = Point2D(cos(s.heel)*getx(p) - sin(s.heel)*gety(p),  sin(s.heel)*getx(p) + cos(s.heel)*gety(p) + s.vshift)
inverse_trans(p, s::BoxShip) = Point2D(cos(s.heel)*getx(p) + sin(s.heel)*(gety(p) - s.vshift), -sin(s.heel)*getx(p) + cos(s.heel)*(gety(p) - s.vshift))

"""
Convert position p to ship coordinates
"""
toship(p, s) = inverse_trans(p, s) - Point2D(0.0u"m", s.height/2 - s.draft)

function corners(s::BoxShip)
  x = s.width/2
  ymax = s.height-s.draft
  ymin = -s.draft
  return forward_trans.(Point2D.([-x,x,x,-x], [ymin,ymin,ymax,ymax]), s)
end

function waterline(s::BoxShip)
  x = 1.4*s.width/2
  y = 0.0u"m"
  return (Point2D(-x,y), Point2D(x,y))
end

function gravitycenter(s::BoxShip)
  return forward_trans(Point2D(0.0u"m", s.KG-s.draft), s)
end

bbox(s::BoxShip) = 1.4 .* (-s.width/2, s.width/2, -s.draft, s.height-s.draft)

function wl_intersect(p1, p2)
  x1 = getx(p1)
  y1 = gety(p1)
  x2 = getx(p2)
  y2 = gety(p2)
  r = (x2-x1)/(y2-y1)
  return Point2D(-y1*r+x1, 0.0u"m")
end

"""
Corners of the underwater part
"""
function carene(s::BoxShip)
  ship_pts = corners(s)
  N = length(ship_pts)
  pts_above = collect(Iterators.filter((p) -> gety(p[2]) > 0u"m", enumerate(ship_pts)))
  sort!(pts_above, lt = (a,b) -> getx(a[2]) < getx(b[2]))
  l1 = pts_above[1][1]
  l2 = l1 % N + 1
  wl_left = wl_intersect(ship_pts[l1], ship_pts[l2])
  r1 = pts_above[end][1]
  r2 = (r1 - 2 + N) % N + 1
  wl_right = wl_intersect(ship_pts[r1], ship_pts[r2])
  pts_below = collect(Iterators.filter((p) -> gety(p[2]) < 0u"m", enumerate(ship_pts)))
  sort!(pts_below, lt = (a,b) -> getx(a[2]) < getx(b[2]))
  return (getindex.(pts_below, 2)..., wl_right, wl_left)
end

function metacenter(s::BoxShip)
  trap = carene(s)
  w_wl = getx(trap[end-1]) - getx(trap[end])
  F = inverse_trans(Point2D((getx(trap[end-1]) + getx(trap[end]))/2, 0.0u"m"),s)
  BM = (w_wl^3 / 12)/carene_area(s)
  B = inverse_trans(buoyancycenter(s),s)
  nvec = F - B
  Mship = B + Point2D(BM*sin(s.heel), BM*cos(s.heel))
  return forward_trans(Mship, s)
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

function triangulate(pts::Union{AbstractArray{ET},NTuple{N,ET} where N}) where ET
  c = centroid(pts)
  triags = Array{NTuple{3,ET}}(length(pts))
  N = length(pts)
  for i in linearindices(pts)
    triags[i] = (c, pts[i], pts[i%N+1])
  end
  return triags
end

function buoyancycenter(s::BoxShip)
  trap = carene(s)
  triags = triangulate(trap)
  areas = area.(triags)
  centroids = centroid.(triags)
  A = sum(areas)
  x = sum(getx.(centroids) .* areas) / A
  y = sum(gety.(centroids) .* areas) / A
  return Point2D(x,y)
end

function carene_area(s::BoxShip)
  trap = carene(s)
  triags = triangulate(trap)
  return sum(area.(triags))
end

function isocarene(s0, θ)
  θ0 = s0.heel
  T = s0.draft
  DT_l = -T
  DT_u = T
  A0 = carene_area(s0)
  for i in 1:50
    DT = (DT_u+DT_l)/2
    s1 = BoxShip(s0.width, s0.height, T, s0.KG, θ+θ0, DT)
    A1 = carene_area(s1)
    if abs(A1-A0)/A0 < 1e-8
      return s1
    end

    if A1 < A0
      DT_u = DT
    else
      DT_l = DT
    end
  end
  error("isocarene did not converge")
end

function labeled_point(p, label, hue="black", radius=5.0)
  sethue(hue)
  circle(p, radius, :fill)
  text(label, Point(p.x+1.2*radius,p.y+1.2*radius), halign=:center, valign=:top)
end

function draw(m::CoordMapping, s::BoxShip, transformation=((p,::BoxShip) -> p); showM=true, showB=true;)
  ship_pts = transformation.(corners(s),s)
  p = remap.(ship_pts, m)
  sethue("black")
  line(p[1], p[2], :stroke)
  line(p[2], p[3], :stroke)
  line(p[4], p[1], :stroke)
  sethue("darkgrey")
  line(p[3], p[4], :stroke)
  labeled_point(remap(transformation(gravitycenter(s),s),m), "G")
  if showM || showB
    M = remap(transformation(metacenter(s),s),m)
    B = remap(transformation(buoyancycenter(s),s),m)
    if showM
      labeled_point(M, "M")
    end
    if showB
      labeled_point(B, "B", "red")
    end
    if showM && showB
      setdash("dot")
      line(M, B, :stroke)
      setdash("solid")
    end
  end
  wl = remap.(transformation.(waterline(s),s), m)
  sethue("blue")
  line(wl[1], wl[2], :stroke)
  #trap = carene(s)
  #F = remap(transformation(Point2D((getx(trap[end-1]) + getx(trap[end]))/2, 0.0u"m"),s),m)
  #labeled_point(F, "F", "blue")
  return m
end

function draw(s::BoxShip, figwidth=300, transformation=((p,::BoxShip) -> p); showM=true, showB=true;)
  draw(drawing(bbox(s), figwidth), s, transformation, showM=showM, showB=showB)
end

setopacity = Luxor.setopacity

end # module
