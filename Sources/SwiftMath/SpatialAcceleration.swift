//
//  SpatialAcceleration.swift
//
//
//  Created by Jacob W Esselstyn on 6/25/23.
//

import Foundation
import Spatial

// MARK: Spatial
public extension RotationAxis3D {
	var tuple: (x: Double, y: Double, z: Double) {
		(x: self.x, y: self.y, z: self.z)
	}
	var tupleCGf: (x: CGFloat, y: CGFloat, z: CGFloat) {
		(x: CGFloat(self.x), y: CGFloat(self.y), z: CGFloat(self.z))
	}
}

import SwiftUI
public extension Angle2D {
	init(_ angle: Angle) {
		self.init(radians: angle.radians)
	}
}
public extension Angle {
	init(_ angle: Angle2D) {
		self.init(radians: angle.radians)
	}
}


// MARK: Accelerate
public extension simd_double3 {
	/// Converts from `Double` to `Float`
	var f: simd_float3 { SIMD3<Float>(self) }
}
extension simd_float3x3: Codable {
	public func encode(to encoder: Encoder) throws {
		var c = encoder.singleValueContainer()
		try c.encode(self)
	}
	public init(from decoder: Decoder) throws {
		let c = try decoder.singleValueContainer()
		self = try c.decode(Self.self)
	}
}
extension simd_double3x3: Codable {
	public func encode(to encoder: Encoder) throws {
		var c = encoder.singleValueContainer()
		try c.encode(self)
	}
	public init(from decoder: Decoder) throws {
		let c = try decoder.singleValueContainer()
		self = try c.decode(Self.self)
	}
}
