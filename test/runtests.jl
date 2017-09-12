using Archimedes
using Base.Test

using Unitful: m, kg

@testset "PointMass" begin
  pm = PointMass(1m, 2m, 3kg, 0.5m)
  @test mass(pm) == 3kg
  @test getx(pm) == 1m
  @test gety(pm) == 2m
  @test radius(pm) == 0.5m
  
  const N = 4
  masses = PointMass.(linspace(1.0m, 4.0m, N), fill(0m, N), fill(1kg, N), fill(1m, N))
  @test bbox(masses) == (1.0m, 4.0m, 0.0m, 0.0m)
  @test centroid(masses) == Point2D(2.5m, 0.0m)
  @test gravitycenter(masses) == PointMass(2.5m, 0.0m, 4.0kg, 2.0m)

  const image_width = 300
  mapping = Archimedes.CoordMapping(image_width, bbox(masses))
  luxor_pt = Archimedes.remap(masses[1], mapping)
  @test luxor_pt.x == -150
  @test luxor_pt.y == 0
end

@testset "BoxShip" begin
  bs = BoxShip(5m, 3m, 2m, 1.5m, 0)
  sc = corners(bs)
  @test sc[1] == Point2D(-2.5m, -2.0m)
  @test sc[3] == Point2D(2.5m, 1.0m)

  triag = Point2D.([0.0m, 1.0m, 0.0m], [0.0m, 0.0m, 2.0m])
  @test Archimedes.area(triag) == 1.0m^2
end