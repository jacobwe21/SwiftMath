//
//  Matrix.swift
//  MySwift
//
//  Created by Jacob W Esselstyn on 1/19/25.
//
import Foundation
import simd
import Accelerate

public struct Matrix<T: BinaryFloatingPoint>: Hashable, CustomStringConvertible {
	
	public let rows: Int
	public let columns: Int
	public private(set) var values: [[T]]
	
	/// Creates a matrix where the outer array is for rows and the inner array is for columns. Throws errors.
	public init(_ values: [[T]], mustBeSquare: Bool) throws {
		guard !values.isEmpty else {
			throw MatrixError.emptyMatrix
		}
		// Verify that the matrix is square
		if mustBeSquare {
			for row in values {
				guard row.count == values.count else {
					throw MatrixError.notSquareMatrix
				}
			}
		}
		self.rows = values.count
		self.columns = values.first!.count
		self.values = values
	}
	/// Creates a matrix where the outer array is for rows and the inner array is for columns.
	public init(_ values: [[T]]) {
		guard !values.isEmpty else {
			fatalError("Matrix Initialization Failure. Matrix is Empty.")
		}
		self.rows = values.count
		self.columns = values.first!.count
		self.values = values
	}
	/// Initialize a matrix of size n x m with a default value
	public init(rows n: Int, columns m: Int, defaultValue: T = 0) {
		self.values = Array(repeating: Array(repeating: defaultValue, count: m), count: n)
		self.rows = n
		self.columns = m
	}
	/// Initialize a matrix of size n x n with a default value
	public init(size n: Int, defaultValue: T = 0) {
		self.values = Array(repeating: Array(repeating: defaultValue, count: n), count: n)
		self.rows = n; self.columns = n
	}
	/// Return the identity matrix of the same size
	public static func identity(size n: Int) -> Matrix<T> {
		var result = Matrix<T>(size: n)
		for i in 0..<n {
			result[i, i] = 1 as! T
		}
		return result
	}
	public init(_ simdMatrix: simd_double2x2) where T == Double {
		self.rows = 2
		self.columns = 2
		let m = simdMatrix.transpose
		self.values = [
			[m.columns.0.x, m.columns.0.y],
			[m.columns.1.x, m.columns.1.y],
		]
	}
	public init(_ simdMatrix: simd_double3x3) where T == Double {
		self.rows = 3
		self.columns = 3
		let m = simdMatrix.transpose
		self.values = [
			[m.columns.0.x, m.columns.0.y, m.columns.0.z],
			[m.columns.1.x, m.columns.1.y, m.columns.1.z],
			[m.columns.2.x, m.columns.2.y, m.columns.2.z],
		]
	}
	public init(_ simdMatrix: simd_double4x4) where T == Double {
		self.rows = 4
		self.columns = 4
		let m = simdMatrix.transpose
		self.values = [
			[m.columns.0.x, m.columns.0.y, m.columns.0.z, m.columns.0.w],
			[m.columns.1.x, m.columns.1.y, m.columns.1.z, m.columns.1.w],
			[m.columns.2.x, m.columns.2.y, m.columns.2.z, m.columns.2.w],
			[m.columns.3.x, m.columns.3.y, m.columns.3.z, m.columns.3.w],
		]
	}
	
	/// Check if rows and columns of each matrix match.
	public static func areDimensionsEqual(lhs: Self, rhs: Self) -> Bool {
		lhs.rows == rhs.rows && lhs.columns == rhs.columns
	}
	
	/// Get a row as an array
	public func row(_ index: Int) -> [T] {
		guard index >= 0 && index < rows else {
			fatalError("Row index out of bounds: \(index)")
		}
		return values[index]
	}
	/// Get a column as an array
	public func column(_ index: Int) -> [T] {
		guard index >= 0 && index < columns else {
			fatalError("Column index out of bounds: \(index)")
		}
		return (0..<rows).map { values[$0][index] }
	}
	
	/// Elementwise addition
	public static func + (lhs: Matrix<T>, rhs: Matrix<T>) throws -> Matrix<T> {
		guard Matrix.areDimensionsEqual(lhs: lhs, rhs: rhs) else { throw MatrixError.nonMatchingDimensions }
		var result = lhs
		for i in 0..<result.rows {
			for j in 0..<result.columns {
				result.values[i][j] = lhs.values[i][j] - rhs.values[i][j]
			}
		}
		return result
	}
	/// Elementwise subtraction
	public static func - (lhs: Matrix<T>, rhs: Matrix<T>) throws -> Matrix<T> {
		guard Matrix.areDimensionsEqual(lhs: lhs, rhs: rhs) else { throw MatrixError.nonMatchingDimensions }
		var result = lhs
		for i in 0..<result.rows {
			for j in 0..<result.columns {
				result.values[i][j] = lhs.values[i][j] - rhs.values[i][j]
			}
		}
		return result
	}
	
	/// Scalar multiplication
	public static func * (lhs: Matrix<T>, rhs: T) -> Matrix<T> {
		var result = lhs
		for i in 0..<result.rows {
			for j in 0..<result.columns {
				result.values[i][j] = result.values[i][j]*rhs
			}
		}
		return result
	}
	/// Scalar division
	public static func / (lhs: Matrix<T>, rhs: T) -> Matrix<T> {
		var result = lhs
		for i in 0..<result.rows {
			for j in 0..<result.columns {
				result.values[i][j] = result.values[i][j]/rhs
			}
		}
		return result
	}
	/// Scalar multiplication
	public static func * (lhs: T, rhs: Matrix<T>) -> Matrix<T> {
		var result = rhs
		for i in 0..<result.rows {
			for j in 0..<result.columns {
				result.values[i][j] = result.values[i][j]*lhs
			}
		}
		return result
	}
	
	/// Elementwise multiplication of two matricies
	public mutating func multiplyElementwise(by other: Self) {
		for i in 0..<rows {
			for j in 0..<columns {
				values[i][j] = values[i][j]*other.values[i][j]
			}
		}
	}
	/// Elementwise multiplication of two matricies
	public func multiplyingElementwise(by other: Self) -> Self {
		var newMatrix = self
		newMatrix.multiplyElementwise(by: other)
		return newMatrix
	}
	
	/// Perform Matrix Multiplication
	public static func * (lhs: Matrix<T>, rhs: Matrix<T>) throws -> Matrix<T> {
		guard lhs.columns == rhs.rows else { throw MatrixError.nonMatchingDimensions }
		var result = Matrix<T>(rows: lhs.rows, columns: rhs.columns)
		for i in 0..<lhs.rows {
			for j in 0..<rhs.columns {
				for k in 0..<lhs.columns {
					result.values[i][j] += lhs.values[i][k] * rhs.values[k][j]
				}
			}
		}
		return result
	}
	
	/// Returns the matrix as a single array, only if matrix is a single column or a single row.
	public func makeVectorArray() throws -> [T] {
		guard rows == 1 || columns == 1 else { throw MatrixError.nonMatchingDimensions }
		if rows == 1 {
			return values[0]
		} else {
			return self.transpose().values[0]
		}
	}
	
	/// Returns the transpose of the matrix (swapping rows for columns)
	public func transpose() -> Self {
		guard rows > 0 && columns > 0 else { return self }
		var result = Matrix(rows: columns, columns: rows)
		for i in 0..<rows {
			for j in 0..<columns {
				result.values[j][i] = self.values[i][j]
			}
		}
		return result
	}
	
	/// Returns `true` if matrix is symmetric
	public var isSymmetric: Bool {
		if isSquare { return self == self.transpose() } else { return false }
	}
	
	/// Check if the matrix is square.
	private var isSquare: Bool { rows == columns }
	
	/// Determinant of the matrix (only defined for square matrices).
	public func determinant() throws -> T {
		guard isSquare else { throw MatrixError.notSquareMatrix }
		return Matrix.calculateDeterminant(for: values)
	}
	
	public func cofactorMatrix() throws -> Matrix<T> {
		guard isSquare else {
			throw MatrixError.notSquareMatrix
		}
		let n = rows
		guard n > 1 else { throw MatrixError.indexOutOfRange }
		if n == 2 {
			return Matrix([
				[values[1][1], -values[1][0]],
				[-values[0][1], values[0][0]]
			])
		}
		var cofactorMatrix = [[T]](repeating: [T](repeating: 0, count: n), count: n)
		for i in 0..<n {
			for j in 0..<n {
				let minor = minorMatrix(removingRow: i, removingColumn: j)
				let sign: T = ((i + j) % 2 == 0) ? 1 : -1
				cofactorMatrix[i][j] = sign * Matrix.calculateDeterminant(for: minor)
			}
		}
		return Matrix(cofactorMatrix)
	}
	public func adjoint() throws -> Matrix<T> {
		// Transpose the cofactor matrix to get the adjoint
		return try self.cofactorMatrix().transpose()
	}
	private func minorMatrix(removingRow rowToRemove: Int, removingColumn columnToRemove: Int) -> [[T]] {
		values.enumerated().compactMap { i, row in
			guard i != rowToRemove else { return nil }
			let filteredRow = row.enumerated().compactMap { j, value in
				j != columnToRemove ? value : nil
			}
			return filteredRow
		}
	}
	
	public func removing(rows rowsToRemove: [Int], columns columnsToRemove: [Int]) throws -> Matrix<T> {
		for row in rowsToRemove {
			if row < 0 || row >= rows { throw MatrixError.indexOutOfRange }
		}
		for column in columnsToRemove {
			if column < 0 || column >= columns { throw MatrixError.indexOutOfRange }
		}
		let rowSet = Set(rowsToRemove)
		let columnSet = Set(columnsToRemove)
		let filteredValues: [[T]] = values.enumerated().compactMap { (i, row) in
			guard !rowSet.contains(i) else { return nil }
			let filteredRow = row.enumerated().compactMap { (j, value) in
				columnSet.contains(j) ? nil : value
			}
			return filteredRow
		}
		return Matrix(filteredValues)
	}
	
	/// Determinant of the matrix (only defined for square matrices).
	public func inverse() throws -> Matrix<T> {
		guard isSquare else { throw MatrixError.notSquareMatrix }
		
		let det = Matrix.calculateDeterminant(for: values)
		if det == 0 { throw MatrixError.singularMatrix }
		
		let size = rows
		if size == 1 { return Matrix([[(1/values[0][0])]]) }
		if size == 3 {
			let m = values

			let a = m[0][0], b = m[0][1], c = m[0][2]
			let d = m[1][0], e = m[1][1], f = m[1][2]
			let g = m[2][0], h = m[2][1], i = m[2][2]

			let A = e * i - f * h
			let B = -(d * i - f * g)
			let C = d * h - e * g
			let D = -(b * i - c * h)
			let E = a * i - c * g
			let F = -(a * h - b * g)
			let G = b * f - c * e
			let H = -(a * f - c * d)
			let I = a * e - b * d
			
			let det = a * A + b * B + c * C
			guard det != 0 else {
				throw MatrixError.singularMatrix
			}
			let resultValues: [[T]] = [
				[A, D, G],
				[B, E, H],
				[C, F, I]
			]
			return Matrix(resultValues)/det
		}
		if size < 5 {
			return (1/det)*(try! self.adjoint())
		}
		
		var augmented = values

		// Append identity matrix to the right of the original
		for i in 0..<size {
			augmented[i] += (0..<size).map { $0 == i ? 1 : 0 }.map(T.init)
		}

		// Perform Gauss-Jordan elimination
		for i in 0..<size {
			// Find pivot
			var pivotRow = i
			for j in i+1..<size where abs(augmented[j][i]) > abs(augmented[pivotRow][i]) {
				pivotRow = j
			}

			// If pivot is zero, matrix is singular
			if augmented[pivotRow][i] == 0 {
				throw MatrixError.singularMatrix
			}

			// Swap rows if needed
			if i != pivotRow {
				augmented.swapAt(i, pivotRow)
			}

			// Normalize pivot row
			let pivot = augmented[i][i]
			for j in 0..<2 * size {
				augmented[i][j] /= pivot
			}

			// Eliminate other rows
			for k in 0..<size where k != i {
				let factor = augmented[k][i]
				for j in 0..<2 * size {
					augmented[k][j] -= factor * augmented[i][j]
				}
			}
		}

		// Extract right half as the inverse
		let inverseValues = augmented.map { Array($0[size..<(2*size)]) }
		return Matrix(inverseValues)
	}
	
	/// Recursive helper function to calculate the determinant.
	private static func calculateDeterminant(for matrix: [[T]]) -> T {
		let n = matrix.count
	
		if n == 1 { // Base case for 1x1 matrix
			return matrix[0][0]
		}
		if n == 2 { // Base case for 2x2 matrix
			return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]
		}
		if n == 3 { // Base case for 3x3 matrix
			let a = matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1])
			let b = matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0])
			let c = matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0])
			return a - b + c
		}
	
		var det: T = 0
		// Calculate determinant using cofactor expansion - O(n!)
		for i in 0..<n {
			// Create the submatrix by excluding the first row and the ith column
			var subMatrix = [[T]]()
			for j in 1..<matrix.count {
				var row = [T]()
				for k in 0..<matrix.count {
					if k != i {
						row.append(matrix[j][k])
					}
				}
				subMatrix.append(row)
			}
			// Recursive determinant calculation with cofactor expansion
			let sign: T = (i % 2 == 0) ? 1 : -1
			det += sign * matrix[0][i] * calculateDeterminant(for: subMatrix)
		}
		// TO-DO
		// Calculate determinant using row reduction (Gaussian elimination) - O(n^3)
		//return rowEcholonForm?.mainDiagonal?.reduce(1, *) ?? 0
		// Calculate determinant using LU decomposition?
		return det
	}
	
	/// Check if indices are valid
	private func isValidIndex(row: Int, col: Int) -> Bool {
		return row >= 0 && row < rows && col >= 0 && col < columns
	}
	/// Fetches or updates a value in the matrix. Will crash if out-of-range.
	public subscript(i: Int, j: Int) -> T {
		get {
			guard isValidIndex(row: i, col: j) else {
				fatalError("Matrix index out of bounds: (\(i), \(j))")
			}
			return values[i][j]
		}
		set(newValue) {
			guard isValidIndex(row: i, col: j) else {
				fatalError("Matrix index out of bounds: (\(i), \(j))")
			}
			values[i][j] = newValue
		}
	}
//	/// Fetches or updates a value in the matrix. Returns `nil` if out-of-range, and does nothing when setting a value if out-of-range or if new value is `nil`.
//	public subscript(i: Int, j: Int) -> T? {
//		get {
//			guard isValidIndex(row: i, col: j) else { return nil }
//			return values[i][j]
//		}
//		set(newValue) {
//			guard isValidIndex(row: i, col: j) && newValue != nil else { return }
//			if let newValue {
//				values[i][j] = newValue
//			}
//		}
//	}
	
	// Error types for matrix operations
	public enum MatrixError: Error {
		case nonMatchingDimensions
		case singularMatrix
		case emptyMatrix
		case notSquareMatrix
		case indexOutOfRange
	}
	
	/// Display row count x column count
	public var sizeDescription: String {
		"\(rows)x\(columns)"
	}
	
	/// Convert matrix to a string representation
	public var description: String {
		var result = "\(rows)x\(columns)\n"
		for i in 0..<rows {
			result += "["
			for j in 0..<columns {
				result += "\(values[i][j])"
				if j < columns - 1 {
					result += ", "
				}
			}
			result += "]\n"
		}
		return result
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
