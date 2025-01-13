//
//  SpatialAcceleration.swift
//
//
//  Created by Jacob W Esselstyn on 6/25/23.
//

import Foundation
import Spatial
import Accelerate
import simd

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

// MARK: Large Sparse Matrices
public extension SparseMatrix_Double {
	
	enum MatrixUpdateError: Error {
		case couldNotFindData
		case invalidRowIndex
		case invalidColumnIndex
	}
	mutating func updateMatrix(at row: Int, column: Int, with newValue: Double) throws {
		guard row >= 0 else { throw MatrixUpdateError.invalidRowIndex }
		guard column >= 0 else { throw MatrixUpdateError.invalidColumnIndex }
		guard row < Int(self.structure.rowCount) else { throw MatrixUpdateError.invalidRowIndex }
		guard column < Int(self.structure.columnCount) else { throw MatrixUpdateError.invalidColumnIndex }
		
		// loop through the matrix to find the index for the data to update
		let nnz = self.structure.columnStarts[Int(self.structure.columnCount)]
		var i = 0
		var currentColumn: Int = 0
		var nextColStarts = self.structure.columnStarts[1]
		while i < nnz {
			if i == nextColStarts {
				currentColumn += 1
				if currentColumn > column { throw MatrixUpdateError.couldNotFindData }
				nextColStarts = self.structure.columnStarts[currentColumn + 1]
			}
			let rowIndex = Int(self.structure.rowIndices[i])
			if rowIndex == row && currentColumn == column {
				self.data[i] = newValue
				return
			}
			i += 1
		}
		throw MatrixUpdateError.couldNotFindData
	}
	
	// https://stackoverflow.com/questions/53006935/how-to-read-values-of-sparsematrix-double-in-swift-4
	var array: [[Double]] {

		let rows = Int(self.structure.rowCount)
		let columns = Int(self.structure.columnCount)
		let nnz = self.structure.columnStarts[Int(self.structure.columnCount)]
				
		var M = Array(repeating: Array(repeating: 0.0, count: columns), count: rows)

		var i = 0
		var currentColumn: Int = 0
		var nextColStarts = self.structure.columnStarts[1]
		while i < nnz {
			if i == nextColStarts {
				currentColumn += 1
				nextColStarts = self.structure.columnStarts[currentColumn + 1]
			}

			let rowIndex = Int(self.structure.rowIndices[i])
			M[rowIndex][currentColumn] = self.data[i]
			i += 1
		}
		return M
	}
	
	subscript (row: Int, column: Int) -> Double {
		get {
			self.array[row][column]
		} set {
			do {
				try self.updateMatrix(at: row, column: column, with: newValue)
			} catch {
				fatalError(error.localizedDescription)
			}
		}
	}
}
//public struct SparseMatrix<T: FloatingPoint> {
//	public var columnCount: Int32
//	public var rowCount: Int32
//	public var columnStarts: [Int]
//	public var rowIndices: [Int32]
//	public var values: [T]
//	public var structure: SparseMatrixStructure
//	
//	public init(_ values: [[T]]) {
//		
//	}
//	
//	public init(rowCount: Int32, columnCount: Int32, columnStarts: [Int], rowIndices: [Int32], values: [T], symmetric: Bool, upperTriangle: Bool = true) {
//		self.columnCount = columnCount
//		self.rowCount = rowCount
//		self.columnStarts = columnStarts
//		self.rowIndices = rowIndices
//		self.values = values
//		
//		var attributes = SparseAttributes_t()
//		if symmetric {
//			if upperTriangle {
//				attributes.triangle = SparseUpperTriangle
//			} else {
//				attributes.triangle = SparseLowerTriangle
//			}
//			attributes.kind = SparseSymmetric
//		}
//		self.columnStarts = columnStarts
//		self.rowIndices = rowIndices
//		structure = self.rowIndices.withUnsafeMutableBufferPointer { rowIndicesPtr in
//			self.columnStarts.withUnsafeMutableBufferPointer { columnStartsPtr in
//				return SparseMatrixStructure(
//					rowCount: rowCount,
//					columnCount: columnCount,
//					columnStarts: columnStartsPtr.baseAddress!,
//					rowIndices: rowIndicesPtr.baseAddress!,
//					attributes: attributes,
//					blockSize: 1
//				)
//			}
//		}
//	}
//}

// MARK: simd Extenions
public extension simd_double2 {
	var magnitude: Double {
		simd_length(self)
	}
}
public extension simd_double3 {
	var magnitude: Double {
		simd_length(self)
	}
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
