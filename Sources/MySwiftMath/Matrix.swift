//
//  Matrix.swift
//  MySwift
//
//  Created by Jacob W Esselstyn on 1/19/25.
//
import Foundation
import simd
import Accelerate
//import Surge

// Error types for matrix operations
public enum MatrixError: String, Error {
	case nonMatchingDimensions = "Dimensions of matrices do not match."
	case singularMatrix = "Matrix is singular."
	case emptyMatrix = "Matrix is empty."
	case notSquareMatrix = "Matrix is not square, and must be square for this operation."
	case indexOutOfRange = "Index is out-of-range."
	case computationFailure = "Computational Failure"
	case matrixNotDecomposable = "Matrix not Decomposable"
}

nonisolated public struct Matrix<T: BinaryFloatingPoint>: Hashable, CustomStringConvertible {
	
	public let rows: Int
	public let columns: Int
	public private(set) var values: [[T]]
	public var flatValues: [T] { values.flatMap(\.self) }
	
	public init<V>(_ matrix: Matrix<V>) where V: BinaryFloatingPoint {
		self.rows = matrix.rows
		self.columns = matrix.columns
		self.values = matrix.values.map({ Array<T>($0.map({T($0)}))})
	}
	
//	@available(iOS 26.0.0, macOS 26.0.0, *)
//	fileprivate init<let bufferSize: Int, V>(_ matrix: InlineMatrix<bufferSize,V>) where V: BinaryFloatingPoint {
//		self.init(rows: matrix.rows, columns: matrix.columns, values: Array(matrix.flatValues))
//	}
//	@available(iOS 26.0.0, macOS 26.0.0, *)
//	fileprivate init<let R: Int, let C: Int, let bufferSize: Int, V>(_ matrix: InlineMatrix<R, C, bufferSize, V>) where V: BinaryFloatingPoint {
//		self.init(rows: R, columns: C, values: Array<V>(matrix.flatValues).map({T($0)}))
//	}
//	@available(iOS 26.0.0, macOS 26.0.0, *)
//	fileprivate init<let N: Int, let bufferSize: Int, V>(_ matrix: InlineSquareMatrix<N, bufferSize, V>) where V: BinaryFloatingPoint {
//		self.init(rows: N, columns: N, values: Array<V>(matrix.flatValues).map({T($0)}))
//	}
	
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
	/// Creates a matrix from a single "vector" array, either as a row or as a column.
	public init(_ values: [T], isRow: Bool) {
		guard !values.isEmpty else {
			fatalError("Matrix Initialization Failure. Matrix is Empty.")
		}
		if isRow {
			self.rows = 1
			self.columns = values.count
			self.values = [values]
		} else {
			self.rows = values.count
			self.columns = 1
			self.values = values.map({[$0]})
		}
	}
	
	/// Creates a matrix from a single array
	/// - Parameters:
	///   - rows: The number of rows
	///   - columns: The number of columns
	///   - values: The data in row-major (row, column) formation
	public init(rows: Int, columns: Int, values: [T]) {
		precondition(values.count == rows * columns)
		self.rows = rows
		self.columns = columns
		self.values = {
			var result: [[T]] = Array(repeating: Array(repeating: 0, count: columns), count: rows)
			for r in 0..<rows {
				for c in 0..<columns {
					result[r][c] = values[r * columns + c]
				}
			}
			return result
		}()
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
	
	public init(rows: Int, columns: Int, _ closure: (_ row: Int, _ column: Int) throws -> T) rethrows {
		var grid: [T] = []
		grid.reserveCapacity(rows * columns)

		for row in 0..<rows {
			for column in 0..<columns {
				grid.append(try closure(row, column))
			}
		}

		self.init(rows: rows, columns: columns, values: grid)
	}

	public static func identity(size: Int) -> Matrix<T> {
		return self.diagonal(rows: size, columns: size, defaultValue: 1.0)
	}
	public static func eyeMatrix(rows: Int, columns: Int) -> Matrix<T> {
		return self.diagonal(rows: rows, columns: columns, defaultValue: 1.0)
	}
	public static func diagonal(rows: Int, columns: Int, defaultValue: T) -> Matrix<T> {
		let count = Swift.min(rows, columns)
		let scalars = repeatElement(defaultValue, count: count)
		return self.diagonal(rows: rows, columns: columns, scalars: scalars)
	}
	public static func diagonal<C>(rows: Int, columns: Int, scalars: C) -> Matrix<T> where C: Collection, C.Element == T {
		var matrix = self.init(rows: rows, columns: columns, defaultValue: 0.0)

		let count = Swift.min(rows, columns)
		precondition(scalars.count == count)

		for (i, scalar) in scalars.enumerated() {
			matrix[i, i] = scalar
		}

		return matrix
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
	
	/// Returns `true` if both matrices are approximately equivalent to 10^-4
	public static func ~= (lhs: Matrix<T>, rhs: Matrix<T>) -> Bool {
		guard lhs.sizeDescription == rhs.sizeDescription else { return false }
		let result = try! lhs - rhs
		for r in 0..<lhs.rows {
			for c in 0..<lhs.columns {
				if abs(result.values[r][c]) > 0.0001 {
					return false
				}
			}
		}
		return true
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
//		var result = Matrix<T>(rows: columns, columns: rows)
//		result.flatValues.withUnsafeMutableBufferPointer { pointer in
//			vDSP_mtransD(flatValues, 1, pointer.baseAddress!, 1, vDSP_Length(columns), vDSP_Length(rows))
//		}
		return result
	}
	
	/// Returns `true` if matrix is symmetric
	public var isSymmetric: Bool {
		if isSquare {
			for r in 0..<rows {
				for c in 0..<columns {
					if self.values[r][c] != self.values[c][r] { return false }
				}
			}
			return true
		} else { return false }
	}
	
	/// Check if the matrix is square.
	private var isSquare: Bool { rows == columns }
	
	/// Determinant of the matrix (only defined for square matrices).
	public func determinant() throws -> T {
		guard isSquare else { throw MatrixError.notSquareMatrix }
		return Self.calculateDeterminant(for: values)
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
	
	/// Recursive helper function to calculate the determinant.
	private static func calculateDeterminant(for matrix: [[T]]) -> T {
		let n = matrix.count
	
		if n == 1 { // 1x1 matrix
			return matrix[0][0]
		}
		if n == 2 { // 2x2 matrix
			return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]
		}
		if n == 3 { // 3x3 matrix
			let a = matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1])
			let b = matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0])
			let c = matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0])
			return a - b + c
		}
		if n < 6 {
			// Calculate determinant using cofactor expansion - O(n!)
			var det: T = 0
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
			return det
		} else {
			// Calculate determinant using LU decomposition - O(n^3)
			do {
				let (L, U, _, rowSwaps) = try Matrix.luDecomposition(of: matrix)
				return calculateDeterminant(fromL: L, withU: U, rowSwaps: rowSwaps)
			} catch {
				if let matrixError = error as? MatrixError {
					if matrixError == MatrixError.singularMatrix { return 0 }
				}
				fatalError(error.localizedDescription)
			}
		}
	}
	/// Calculate the determinant from LU Decomposition
	private static func calculateDeterminant(fromL L: [[T]], withU U: [[T]], rowSwaps: Int) -> T {
		//let detL: T = 1
		var detU: T = 1
		// Product of diagonal elements
		for i in 0..<L.count {
			//detL *= L[i][i] // Diagonal elements are all 1.
			detU *= U[i][i]
		}
		//let det = detL*detU
		// Adjust sign for number of row swaps
		return (rowSwaps % 2 == 0) ? detU : -detU
	}
	
	public static func solve(A: Matrix<T>, b: Matrix<T>) throws -> Matrix<T> {
		guard A.isSquare else { throw MatrixError.notSquareMatrix }
		guard b.columns == 1 && b.rows == A.rows else { throw MatrixError.nonMatchingDimensions }
		let n = A.rows
		if n < 6 {
			let inverseA = try A.inverse()
			return try inverseA * b
		}
		
		let (L, U, P, _) = try Matrix.luDecomposition(of: A.values)
		
		let x = try Matrix.solve(L: L, U: U, P: P, b: b.flatValues)
		return Matrix(x, isRow: false)
	}
	public static func solve(L: [[T]], U: [[T]], P: [Int], b: [T]) throws -> [T] {
		guard L.count == L.first?.count && U.count == U.first?.count else {
			throw MatrixError.notSquareMatrix
		}
		guard b.count == L.count && L.count == U.count && P.count == L.count else {
			throw MatrixError.nonMatchingDimensions
		}
		let n = L.count
		
		// Helper to apply permutation to identity column
		func permute(_ e: [T]) -> [T] {
			var result = Array<T>(repeating: 0.0, count: n)
			for i in 0..<n {
				result[i] = e[P[i]]
			}
			return result
		}
		
		let Pb = permute(b)
		
		// Forward substitution: solve L * y = Pb
		var y = Array<T>(repeating: 0.0, count: n)
		for j in 0..<n {
			y[j] = Pb[j]
			for k in 0..<j {
				y[j] -= L[j][k] * y[k]
			}
			// L diagonal is 1, so divide unnecessary
		}
		
		// Backward substitution: solve U * x = y
		var x = Array<T>(repeating: 0.0, count: n)
		for j in stride(from: n - 1, through: 0, by: -1) {
			x[j] = y[j]
			for k in (j + 1)..<n {
				x[j] -= U[j][k] * x[k]
			}
			guard U[j][j] != 0 else { throw MatrixError.singularMatrix } // Singular matrix
			x[j] /= U[j][j]
		}
		return x
	}
	
	/// Determinant of the matrix (only defined for square matrices).
	public func inverse() throws -> Matrix<T> {
		guard isSquare else { throw MatrixError.notSquareMatrix }
		let n = rows
		if n < 6 {
			let det = Matrix.calculateDeterminant(for: values)
			if det == 0 { throw MatrixError.singularMatrix }
			
			if n == 1 { return Matrix([[(1/values[0][0])]]) }
			if n == 3 {
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
				
				//let det = a * A + b * B + c * C
				let resultValues: [[T]] = [
					[A, D, G],
					[B, E, H],
					[C, F, I]
				]
				return Matrix(resultValues)/det
			}
			// Use Cofactor Expansion for small matrices (n! < n^3).
			return (1/det)*(try! self.adjoint())
		}
		
		// TODO: Implement Cholesky decomposition for symmetric positive definite matrices
		// Use LU Decomposition
		let (L, U, P, _) = try Matrix.luDecomposition(of: self.values)
		var invCols: [[T]] = []
		
		// Helper to apply permutation to identity column
		func permute(_ e: [T]) -> [T] {
			var result = Array<T>(repeating: 0.0, count: n)
			for i in 0..<n {
				result[i] = e[P[i]]
			}
			return result
		}
		
		// Solve for each column of the inverse
		for i in 0..<n {
			// Unit vector (column of identity)
			var e = Array<T>(repeating: 0.0, count: n)
			e[i] = 1.0
			let b = permute(e)
			
			// Forward substitution: solve L * y = b
			var y = Array<T>(repeating: 0.0, count: n)
			for j in 0..<n {
				y[j] = b[j]
				for k in 0..<j {
					y[j] -= L[j][k] * y[k]
				}
				// L diagonal is 1, so divide unnecessary
			}
			
			// Backward substitution: solve U * x = y
			var x = Array<T>(repeating: 0.0, count: n)
			for j in stride(from: n - 1, through: 0, by: -1) {
				x[j] = y[j]
				for k in (j + 1)..<n {
					x[j] -= U[j][k] * x[k]
				}
				guard U[j][j] != 0 else { throw MatrixError.singularMatrix } // Singular matrix
				x[j] /= U[j][j]
			}
			
			invCols.append(x)
		}
		
		// Transpose columns to rows
		let inverseValues = (0..<n).map { i in
			invCols.map { $0[i] }
		}
		return Matrix(inverseValues)

		// Perform Gauss-Jordan elimination
//		var augmented = values
//
//		// Append identity matrix to the right of the original
//		for i in 0..<size {
//			augmented[i] += (0..<size).map { $0 == i ? 1 : 0 }.map(T.init)
//		}
//
//		for i in 0..<size {
//			// Find pivot
//			var pivotRow = i
//			for j in i+1..<size where abs(augmented[j][i]) > abs(augmented[pivotRow][i]) {
//				pivotRow = j
//			}
//
//			// If pivot is zero, matrix is singular
//			if augmented[pivotRow][i] == 0 {
//				throw MatrixError.singularMatrix
//			}
//
//			// Swap rows if needed
//			if i != pivotRow {
//				augmented.swapAt(i, pivotRow)
//			}
//
//			// Normalize pivot row
//			let pivot = augmented[i][i]
//			for j in 0..<2 * size {
//				augmented[i][j] /= pivot
//			}
//
//			// Eliminate other rows
//			for k in 0..<size where k != i {
//				let factor = augmented[k][i]
//				for j in 0..<2 * size {
//					augmented[k][j] -= factor * augmented[i][j]
//				}
//			}
//		}
//
//		// Extract right half as the inverse
//		let inverseValues = augmented.map { Array($0[size..<(2*size)]) }
//		return Matrix(inverseValues)
	}
	
	/// Returns (L, U, P, swapCount) for LU decomposition with partial pivoting.
	/// - Returns: Optional tuple (L, U, P, swapCount), or nil if decomposition fails.
	public static func luDecomposition(of matrix: [[T]], tolerance: T = 0.000001) throws -> (L: [[T]], U: [[T]], P: [Int], swapCount: Int) {
		guard matrix.count == matrix.first?.count else { throw MatrixError.notSquareMatrix }
		let n = matrix.count
		var A = matrix
		var L = Array(repeating: Array<T>(repeating: 0.0, count: n), count: n)
		var U = Array(repeating: Array<T>(repeating: 0.0, count: n), count: n)
		var P: Array<Int> = Array(0..<n)
		var swapCount = 0
		
		let swapsAreRequired: Bool = {
			for i in 0..<n {
				if abs(A[i][i]) < tolerance { return true }
			}
			return false
		}()
		
		if swapsAreRequired {
			for i in 0..<n {
				// Partial Pivoting
				var maxRow = i
				var maxVal: T = 0.0
				for k in i..<n {
					let kABS = abs(A[k][i])
					if kABS > maxVal {
						maxVal = kABS
						maxRow = k
					}
				}
				
				// Zero pivot => singular (or near singular) matrix
				if maxVal < tolerance {
					throw MatrixError.singularMatrix
				}
				
				// Swap rows in A and P
				if maxRow != i {
					A.swapAt(i, maxRow) // Pivot A
					P.swapAt(i, maxRow) // Pivot P
					swapCount += 1
				}
			}
		}

		// Doolittle algorithm + Partial Pivoting
		for i in 0..<n {
			
			// Upper Triangular
			for j in i..<n {
				// Summation of L(i, k) * U(k, j)
				var sum: T = 0.0
				for k in 0..<i {
					sum += L[i][k] * U[k][j]
				}
				// Evaluating U(i, j)
				U[i][j] = A[i][j] - sum
			}

			
			// Lower Triangular
			L[i][i] = 1.0 // unit diagonal
			for j in (i+1)..<n {
				// Summation of L(j, k) * U(k, i)
				var sum: T = 0.0
				for k in 0..<i {
					sum += L[j][k] * U[k][i]
				}
				// Evaluating L(j, i)
				if U[i][i] == 0.0 { throw MatrixError.singularMatrix } // Singular matrix
				L[j][i] = (A[j][i] - sum) / U[i][i]
			}
		}
		return (L, U, P, swapCount)
	}
	/// Builds a permutation matrix from a row permutation array.
	/// Each index `i` in `permutation` maps original row `i` to `permutation[i]`.
	public static func permutationMatrix(from permutation: [Int]) -> Matrix<T> {
		let n = permutation.count
		var values = Array(repeating: Array<T>(repeating: 0.0, count: n), count: n)
		for i in 0..<n {
			values[permutation[i]][i] = 1.0
		}
		return Matrix(values)
	}
	
	// UNTESTED - written by AI - LU decomposition is about 2x faster than QR usually for solving matrices, but QR is best for eigenvalues and eigenvectors.
	private func qrDecompose(iterations: Int = 10) -> (Q: [[T]], R: [[T]]) {
		var A = self.values
		let n = A.count
		let m = A[0].count
		var Q: [[T]] = Array(repeating: Array(repeating: 0.0, count: m), count: n)
		var R: [[T]] = Array(repeating: Array(repeating: 0.0, count: m), count: m)

		for k in 0..<m {
			var norm: T = 0.0
			for i in 0..<n {
				norm += A[i][k] * A[i][k]
			}
			norm = sqrt(norm)
			for i in 0..<n {
				Q[i][k] = A[i][k] / norm
			}
			R[k][k] = norm

			for j in (k+1)..<m {
				var dot: T = 0.0
				for i in 0..<n {
					dot += Q[i][k] * A[i][j]
				}
				R[k][j] = dot
				for i in 0..<n {
					A[i][j] -= Q[i][k] * dot
				}
			}
		}
		return (Q, R)
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
	
	// MARK: Surge
	// Surge MIT License:

	//Copyright © 2014-2019 the Surge contributors
	//
	//Permission is hereby granted, free of charge, to any person obtaining a copy
	//of this software and associated documentation files (the "Software"), to deal
	//in the Software without restriction, including without limitation the rights
	//to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
	//copies of the Software, and to permit persons to whom the Software is
	//furnished to do so, subject to the following conditions:
	//
	//The above copyright notice and this permission notice shall be included in
	//all copies or substantial portions of the Software.
	//
	//THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
	//IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
	//FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
	//AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
	//LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
	//OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
	//THE SOFTWARE.
	
	// ALTERNATE Subscripts for flatValues formulation.
//	public subscript(row: Int, column: Int) -> T {
//		get {
//			assert(isValidIndex(row: row, col: column))
//			return flatValues[(row * columns) + column]
//		}
//		set {
//			assert(isValidIndex(row: row, col: column))
//			flatValues[(row * columns) + column] = newValue
//		}
//	}
//	public subscript(row row: Int) -> [T] {
//		get {
//			assert(row < rows)
//			let startIndex = row * columns
//			let endIndex = row * columns + columns
//			return Array(grid[startIndex..<endIndex])
//		}
//
//		set {
//			assert(row < rows)
//			assert(newValue.count == columns)
//			let startIndex = row * columns
//			let endIndex = row * columns + columns
//			grid.replaceSubrange(startIndex..<endIndex, with: newValue)
//		}
//	}
//	public subscript(column column: Int) -> [T] {
//		get {
//			var result = [T](repeating: 0.0, count: rows)
//			for i in 0..<rows {
//				let index = i * columns + column
//				result[i] = self.grid[index]
//			}
//			return result
//		}
//
//		set {
//			assert(column < columns)
//			assert(newValue.count == rows)
//			for i in 0..<rows {
//				let index = i * columns + column
//				grid[index] = newValue[i]
//			}
//		}
//	}
	
	public func eigenvalues(maxIterations: Int = 100, tolerance: Double = 1e-10) throws -> [(x: T, i: T)] {
		guard isSquare else { throw MatrixError.notSquareMatrix }
		if rows == 1 {
			return values[0][0].isNaN ? [] : [((values[0][0]), T(0))]
		}
		if rows == 2 {
			let a = values[0][0]
			let b = values[0][1]
			let c = values[1][0]
			let d = values[1][1]

			let trace = a + d
			let determinant = a * d - b * c
			let discriminant = trace * trace - 4 * determinant

			if discriminant >= 0 {
				let sqrtDiscriminant = discriminant.squareRoot()
				let lambda1 = (trace + sqrtDiscriminant) / 2
				let lambda2 = (trace - sqrtDiscriminant) / 2
				return [(lambda1,0.0), (lambda2,0.0)]
			}
		}
		return try eigenDecompose(computeEigenVectors: false).eigenValues
	}
	
	/// Decomposes a square matrix into its eigenvalues and left and right eigenvectors.
	/// The decomposition may result in complex numbers, represented by (T, T), which
	///   are the (real, imaginary) parts of the complex number.
	/// - Returns: a struct with the eigen values and left and right eigen vectors using (T, T)
	///   to represent a complex number.
	func eigenDecompose(computeEigenVectors: Bool) throws -> MatrixEigenDecompositionResult<T> {
		let input = Matrix<Double>(self)
		let decomposition = try eigenDecomposeOLD(input, computeEigenVectors: computeEigenVectors)
		
//		guard rows == columns else {
//			throw MatrixError.notSquareMatrix
//		}
//		let (Q,R) = self.qrDecompose()
//		let lambdas: [(x: T, i: T)] = []
//		var rightEigenVectors: [[(T, T)]] = [[]]
//		if computeEigenVectors  {
//			for lambda in lambdas {
//				var A = try self-Matrix<T>.identity(size: rows)*lambda.x
//				let realVector = try Matrix.solve(A: A, b: Matrix(Array<T>(repeating: 0, count: rows), isRow: false))
//				A = try self-Matrix<T>.identity(size: rows)*lambda.i
//				let imaginaryVector = try Matrix.solve(A: A, b: Matrix(Array<T>(repeating: 0, count: rows), isRow: false))
//				var eigenVector: [(T, T)] = []
//				for i in 0..<rows {
//					eigenVector.append((realVector.flatValues[i],imaginaryVector.flatValues[i]))
//				}
//				rightEigenVectors.append(eigenVector)
//			}
//		}
		
		return MatrixEigenDecompositionResult<T>(
			eigenValues: decomposition.eigenValues.map { (T($0.0), T($0.1)) },
			leftEigenVectors: decomposition.leftEigenVectors.map { $0.map { (T($0.0), T($0.1)) } },
			rightEigenVectors: decomposition.rightEigenVectors.map { $0.map { (T($0.0), T($0.1)) } }
		)
	}

	/// Decomposes a square matrix into its eigenvalues and left and right eigenvectors.
	/// The decomposition may result in complex numbers, represented by (Double, Double), which
	///   are the (real, imaginary) parts of the complex number.
	/// - Parameters:
	///   - lhs: a square matrix
	/// - Returns: a struct with the eigen values and left and right eigen vectors using (Double, Double)
	///   to represent a complex number.
	func eigenDecomposeOLD(_ lhs: Matrix<Double>, computeEigenVectors: Bool) throws -> MatrixEigenDecompositionResult<Double> {
		
		guard lhs.rows == lhs.columns else {
			throw MatrixError.notSquareMatrix
		}

		// dgeev_ needs column-major matrices, so transpose 'lhs'.
		var matrixGrid: [__CLPK_doublereal] = lhs.transpose().flatValues
		var matrixRowCount = __CLPK_integer(lhs.rows)
		let matrixColCount = matrixRowCount
		var eigenValueCount = matrixRowCount
		var leftEigenVectorCount = matrixRowCount
		var rightEigenVectorCount = matrixRowCount

		var workspaceQuery: Double = 0.0
		var workspaceSize = __CLPK_integer(-1)
		var error: __CLPK_integer = 0

		var eigenValueRealParts = [Double](repeating: 0, count: Int(eigenValueCount))
		var eigenValueImaginaryParts = [Double](repeating: 0, count: Int(eigenValueCount))
		var leftEigenVectorWork = [Double](repeating: 0, count: Int(leftEigenVectorCount * matrixColCount))
		var rightEigenVectorWork = [Double](repeating: 0, count: Int(rightEigenVectorCount * matrixColCount))

		var decompositionJobVL: [CChar]
		var decompositionJobVR: [CChar]
		
		if computeEigenVectors  {
			decompositionJobVL = [0x56, 0x00] // "V" (compute)
			decompositionJobVR = [0x56, 0x00] // "V" (compute)
		} else {
			decompositionJobVL = Array("N".utf8CString) // "N" (do not compute)
			decompositionJobVR = Array("N".utf8CString) // "N" (do not compute)
		}

		// Call dgeev to find out how much workspace to allocate
		dgeev_(&decompositionJobVL, &decompositionJobVR, &matrixRowCount, &matrixGrid, &eigenValueCount, &eigenValueRealParts, &eigenValueImaginaryParts, &leftEigenVectorWork, &leftEigenVectorCount, &rightEigenVectorWork, &rightEigenVectorCount, &workspaceQuery, &workspaceSize, &error)
		if error != 0 {
			throw MatrixError.matrixNotDecomposable
		}

		// Allocate the workspace and call dgeev again to do the actual decomposition
		var workspace = [Double](repeating: 0.0, count: Int(workspaceQuery))
		workspaceSize = __CLPK_integer(workspaceQuery)
		dgeev_(&decompositionJobVL, &decompositionJobVR, &matrixRowCount, &matrixGrid, &eigenValueCount, &eigenValueRealParts, &eigenValueImaginaryParts, &leftEigenVectorWork, &leftEigenVectorCount, &rightEigenVectorWork, &rightEigenVectorCount, &workspace, &workspaceSize, &error)
		if error != 0 {
			throw MatrixError.matrixNotDecomposable
		}

		return MatrixEigenDecompositionResult<Double>(rowCount: lhs.rows, eigenValueRealParts: eigenValueRealParts, eigenValueImaginaryParts: eigenValueImaginaryParts, leftEigenVectorWork: leftEigenVectorWork, rightEigenVectorWork: rightEigenVectorWork)
	}
	
	/// Holds the result of eigendecomposition. The (Scalar, Scalar) used
	/// in the property types represents a complex number with (real, imaginary) parts.
	struct MatrixEigenDecompositionResult<Scalar: BinaryFloatingPoint> {
		public let eigenValues: [(Scalar, Scalar)]
		public let leftEigenVectors: [[(Scalar, Scalar)]]
		public let rightEigenVectors: [[(Scalar, Scalar)]]
		
		public init(eigenValues: [(Scalar, Scalar)], leftEigenVectors: [[(Scalar, Scalar)]], rightEigenVectors: [[(Scalar, Scalar)]]) {
			self.eigenValues = eigenValues
			self.leftEigenVectors = leftEigenVectors
			self.rightEigenVectors = rightEigenVectors
		}
		
		public init(rowCount: Int, eigenValueRealParts: [Scalar], eigenValueImaginaryParts: [Scalar], leftEigenVectorWork: [Scalar], rightEigenVectorWork: [Scalar]) {
			// The eigenvalues are an array of (real, imaginary) results from dgeev
			self.eigenValues = Array(zip(eigenValueRealParts, eigenValueImaginaryParts))
			
			// Build the left and right eigenvectors
			let emptyVector = [(Scalar, Scalar)](repeating: (0.0, 0.0), count: rowCount)
			var leftEigenVectors = [[(Scalar, Scalar)]](repeating: emptyVector, count: rowCount)
			MatrixEigenDecompositionResult.buildEigenVector(eigenValueImaginaryParts: eigenValueImaginaryParts, eigenVectorWork: leftEigenVectorWork, result: &leftEigenVectors)
			
			var rightEigenVectors = [[(Scalar, Scalar)]](repeating: emptyVector, count: rowCount)
			MatrixEigenDecompositionResult.buildEigenVector(eigenValueImaginaryParts: eigenValueImaginaryParts, eigenVectorWork: rightEigenVectorWork, result: &rightEigenVectors)
			
			self.leftEigenVectors = leftEigenVectors
			self.rightEigenVectors = rightEigenVectors
		}
		
		// Convert the result of dgeev into an array of complex numbers
		// See Intel's documentation on column-major results for sample code that this is based on:
		// https://software.intel.com/sites/products/documentation/doclib/mkl_sa/11/mkl_lapack_examples/dgeev.htm
		private static func buildEigenVector<S>(eigenValueImaginaryParts: [S], eigenVectorWork: [S], result: inout [[(S, S)]]) where S: FloatingPoint & ExpressibleByFloatLiteral {
			// row and col count are the same because result must be square.
			let rowColCount = result.count

			for row in 0..<rowColCount {
				var col = 0
				while col < rowColCount {
					if eigenValueImaginaryParts[col] == 0.0 {
						// v is column-major
						result[row][col] = (eigenVectorWork[row + rowColCount * col], 0.0)
						col += 1
					} else {
						// v is column-major
						result[row][col] = (eigenVectorWork[row + col * rowColCount], eigenVectorWork[row + rowColCount * (col + 1)])
						result[row][col + 1] = (eigenVectorWork[row + col * rowColCount], -eigenVectorWork[row + rowColCount * (col + 1)])
						col += 2
					}
				}
			}
		}

	}
}
//extension Matrix: Collection {
//	public subscript(_ row: Int) -> ArraySlice<T> {
//		let startIndex = row * columns
//		let endIndex = startIndex + columns
//		return self.flatValues[startIndex..<endIndex]
//	}
//	public var startIndex: Int { return 0 }
//	public var endIndex: Int { return self.rows }
//	public func index(after i: Int) -> Int {
//		return i + 1
//	}
//}

//@available(iOS 26.0.0, macOS 26.0.0, *)
//extension Array where Element: BinaryFloatingPoint {
//	init<let bufferSize: Int>(_ inlineArray: borrowing InlineArray<bufferSize, Element>) {
//		var result = Array<Element>()
//		for i in inlineArray.indices {
//			result.append(inlineArray[i])
//		}
//		self = result
////		self = Array(unsafeUninitializedCapacity: bufferSize) { pointerBuffer, initializedCount in
////			initializedCount = bufferSize
////			for i in inlineArray.indices {
////				pointerBuffer[i] = inlineArray[i]
////				//pointerBuffer.baseAddress.?.initialize(to: flatValues[i])
////			}
////		}
//	}
//}

//@available(iOS 26.0.0, macOS 26.0.0, *)
//extension InlineArray: @retroactive Equatable where Element: Equatable {
//	public static func == <let lhsN: Int, let rhsN: Int>(lhs: InlineArray<lhsN, Self.Element>, rhs: InlineArray<rhsN, Self.Element>) -> Bool {
//		guard lhs.count == rhs.count else { return false }
//		for i in lhs.indices {
//			if lhs[i] != rhs[i] { return false }
//		}
//		return true
//	}
//}

//@available(iOS 26.0.0, macOS 26.0.0, *)
//nonisolated private struct InlineMatrix<let buffer: Int, T: BinaryFloatingPoint>: Equatable, CustomStringConvertible, Copyable {
//	
//	public let rows: Int
//	public let columns: Int
//	public var values: [[T]] {
//		get {
//			var result = [[T]]()
//			for i in 0..<buffer {
//				let r = i / columns
//				let c = i % columns
//				result[r][c] = flatValues[i]
//			}
//			return result
//		}
//		set {
//			for i in 0..<buffer {
//				let r = i / columns
//				let c = i % columns
//				flatValues[i] = newValue[r][c]
//			}
//		}
//	}
//	public private(set)var flatValues: InlineArray<buffer, T>
//	
//	static func == (lhs: InlineMatrix<buffer, T>, rhs: InlineMatrix<buffer, T>) -> Bool {
//		lhs.rows == rhs.rows && lhs.columns == rhs.columns && lhs.flatValues == rhs.flatValues
//	}
//	func isIdentical(to rhs: InlineMatrix<buffer, T>) -> Bool {
//		self.rows == rhs.rows && self.columns == rhs.columns && self.flatValues.span.isIdentical(to: rhs.flatValues.span)
//	}
//	
//	public init<V>(_ matrix: InlineMatrix<buffer, V>) where V: BinaryFloatingPoint {
//		self.rows = matrix.rows
//		self.columns = matrix.columns
//		self.flatValues = InlineArray<bufferSize, T>() { T(matrix.flatValues[$0]) }
//	}
//	/// Creates a matrix where the outer array is for rows and the inner array is for columns. Throws errors.
//	public init(_ values: [[T]], mustBeSquare: Bool) throws {
//		guard !values.isEmpty else {
//			throw MatrixError.emptyMatrix
//		}
//		// Verify that the matrix is square
//		if mustBeSquare {
//			for row in values {
//				guard row.count == values.count else {
//					throw MatrixError.notSquareMatrix
//				}
//			}
//		}
//		self.rows = values.count
//		self.columns = values.first!.count
//		let tempFlatValues = values.flatMap(\.self)
//		self.flatValues = InlineArray() {
//			tempFlatValues[$0]
//		}
//	}
//	/// Creates a matrix where the outer array is for rows and the inner array is for columns.
//	public init(_ values: [[T]]) {
//		guard !values.isEmpty else {
//			fatalError("Matrix Initialization Failure. Matrix is Empty.")
//		}
//		self.rows = values.count
//		self.columns = values.first!.count
//		let tempFlatValues = values.flatMap(\.self)
//		self.flatValues = InlineArray() {
//			tempFlatValues[$0]
//		}
//	}
//	/// Creates a matrix from a single "vector" array, either as a row or as a column.
//	public init(_ values: [T], isRow: Bool) {
//		guard !values.isEmpty else {
//			fatalError("Matrix Initialization Failure. Matrix is Empty.")
//		}
//		if isRow {
//			self.rows = 1
//			self.columns = values.count
//			self.flatValues = InlineArray() { values[$0] }
//		} else {
//			self.rows = values.count
//			self.columns = 1
//			self.flatValues = InlineArray() { values[$0] }
//		}
//	}
//	
//	/// Creates a matrix from a single array
//	/// - Parameters:
//	///   - columns: The number of columns
//	///   - values: The data in row-major (row, column) formation
//	public init(columns m: Int, values: [T]) {
//		precondition(buffer%m == 0)
//		self.rows = buffer/m
//		self.columns = m
//		self.flatValues = InlineArray<buffer, T>() { values[$0] }
//	}
//	/// Creates a matrix from a single array
//	/// - Parameters:
//	///   - columns: The number of columns
//	///   - data: The data in row-major (row, column) formation
//	public init(columns m: Int, data: InlineArray<buffer, T>) {
//		precondition(buffer%m == 0)
//		self.rows = buffer/m
//		self.columns = m
//		self.flatValues = data
//	}
//	/// Initialize a matrix of size n x m with a default value
//	public init(columns m: Int, defaultValue: T = 0) {
//		precondition(buffer%m == 0)
//		self.rows = buffer/m
//		self.columns = m
//		self.flatValues = InlineArray<buffer, T>(repeating: defaultValue)
//	}
//	/// Initialize a matrix of size n x n with a default value
//	public init(size n: Int? = nil, defaultValue: T = 0) {
//		self.flatValues = InlineArray<buffer, T>(repeating: defaultValue)
//		self.rows = n ?? Int(sqrt(Double(buffer)))
//		self.columns = n ?? Int(sqrt(Double(buffer)))
//	}
//	
//	public init(columns m: Int, _ closure: (_ row: Int, _ column: Int) throws -> T) rethrows {
//		precondition(buffer%m == 0)
//		self.flatValues = try InlineArray<buffer, T>() { i in
//			let r = i / m
//			let c = i % m
//			return try closure(r, c)
//		}
//		rows = buffer/m
//		columns = m
//	}
//	public init(columns m: Int, _ closure: (_ i: Int) throws -> T) rethrows {
//		precondition(buffer%m == 0)
//		self.flatValues = try InlineArray<buffer, T>() { i in
//			return try closure(i)
//		}
//		rows = buffer/m
//		columns = m
//	}
//
//	public static func identity(size: Int? = nil) -> InlineMatrix<buffer, T> {
//		let n = size ?? Int(sqrt(Double(buffer)))
//		return InlineMatrix<buffer, T>(columns: n) { i in
//			if i / n == i % n { 1 } else { 0 }
//		}
//	}
//	public static func diagonal(size: Int? = nil, defaultValue: T) -> InlineMatrix<buffer, T> {
//		let n = size ?? Int(sqrt(Double(buffer)))
//		return InlineMatrix<buffer, T>(columns: n) { i in
//			if i / n == i % n { defaultValue } else { 0 }
//		}
//	}
//	public static func diagonal<C>(size: Int? = nil, scalars: C) -> InlineMatrix<buffer, T> where C: Collection, C.Element == T {
//		let n = size ?? Int(sqrt(Double(buffer)))
//		var iScalar: C.Index = scalars.index(scalars.startIndex, offsetBy: -1)
//		return InlineMatrix<buffer, T>(columns: n) { i in
//			if i / n == i % n {
//				iScalar = scalars.index(after: iScalar)
//				return scalars[iScalar]
//			} else { return 0 }
//		}
//	}
//
//	public init(_ simdMatrix: simd_double2x2) where T == Double {
//		self.rows = 2
//		self.columns = 2
//		let m = simdMatrix.transpose
//		let temp = [
//			m.columns.0.x, m.columns.0.y,
//			m.columns.1.x, m.columns.1.y,
//		]
//		self.flatValues = InlineArray() { temp[$0] }
//	}
//	public init(_ simdMatrix: simd_double3x3) where T == Double {
//		self.rows = 3
//		self.columns = 3
//		let m = simdMatrix.transpose
//		let temp = [
//			m.columns.0.x, m.columns.0.y, m.columns.0.z,
//			m.columns.1.x, m.columns.1.y, m.columns.1.z,
//			m.columns.2.x, m.columns.2.y, m.columns.2.z,
//		]
//		self.flatValues = InlineArray() { temp[$0] }
//	}
//	public init(_ simdMatrix: simd_double4x4) where T == Double {
//		self.rows = 4
//		self.columns = 4
//		let m = simdMatrix.transpose
//		let temp = [
//			m.columns.0.x, m.columns.0.y, m.columns.0.z, m.columns.0.w,
//			m.columns.1.x, m.columns.1.y, m.columns.1.z, m.columns.1.w,
//			m.columns.2.x, m.columns.2.y, m.columns.2.z, m.columns.2.w,
//			m.columns.3.x, m.columns.3.y, m.columns.3.z, m.columns.3.w,
//		]
//		self.flatValues = InlineArray() { temp[$0] }
//	}
//	
//	/// Check if rows and columns of each matrix match.
//	public static func areDimensionsEqual(lhs: Self, rhs: Self) -> Bool {
//		lhs.rows == rhs.rows && lhs.columns == rhs.columns
//	}
//	
////	public subscript(row row: Int) -> Array<T> {
////		get {
////			assert(row < rows)
////			let startIndex = row * columns
////			let endIndex = row * columns + columns
////			return Array(flatValues[startIndex..<endIndex])
////		}
////		set {
////			assert(row < rows)
////			assert(newValue.count == columns)
////			let startIndex = row * columns
////			let endIndex = row * columns + columns
////			flatValues.replaceSubrange(startIndex..<endIndex, with: newValue)
////		}
////	}
////	public subscript(column column: Int) -> Array<T> {
////		get {
////			var result = [T](repeating: 0.0, count: rows)
////			for i in 0..<rows {
////				let index = i * columns + column
////				result[i] = self.flatValues[index]
////			}
////			return result
////		}
////		set {
////			assert(column < columns)
////			assert(newValue.count == rows)
////			for i in 0..<rows {
////				let index = i * columns + column
////				flatValues[index] = newValue[i]
////			}
////		}
////	}
//	
//	/// Returns `true` if both matrices are approximately equivalent to 10^-4
//	public static func ~= (lhs: InlineMatrix<buffer,T>, rhs: InlineMatrix<buffer,T>) -> Bool {
//		guard InlineMatrix.areDimensionsEqual(lhs: lhs, rhs: rhs) else { return false }
//		let result = try! lhs - rhs
//		for i in 0..<buffer {
//			if abs(result.flatValues[i]) > 0.000001 {
//				return false
//			}
//		}
//		return true
//	}
//	
//	/// Elementwise addition
//	public static func + (lhs: InlineMatrix<buffer,T>, rhs: InlineMatrix<buffer,T>) throws -> InlineMatrix<buffer,T> {
//		guard InlineMatrix.areDimensionsEqual(lhs: lhs, rhs: rhs) else { throw MatrixError.nonMatchingDimensions }
//		var result = lhs
//		for i in 0..<buffer {
//			result.flatValues[i] = lhs.flatValues[i] - rhs.flatValues[i]
//		}
//		return result
//	}
//	/// Elementwise subtraction
//	public static func - (lhs: InlineMatrix<buffer,T>, rhs: InlineMatrix<buffer,T>) throws -> InlineMatrix<buffer,T> {
//		guard InlineMatrix.areDimensionsEqual(lhs: lhs, rhs: rhs) else { throw MatrixError.nonMatchingDimensions }
//		var result = lhs
//		for i in 0..<buffer {
//			result.flatValues[i] = lhs.flatValues[i] - rhs.flatValues[i]
//		}
//		return result
//	}
//	
//	/// Scalar multiplication
//	public static func * (lhs: InlineMatrix<buffer,T>, rhs: T) -> InlineMatrix<buffer,T> {
//		var result = lhs
//		for i in 0..<buffer {
//			result.flatValues[i] = lhs.flatValues[i] * rhs
//		}
//		return result
//	}
//	/// Scalar division
//	public static func / (lhs: InlineMatrix<buffer,T>, rhs: T) -> InlineMatrix<buffer,T> {
//		var result = lhs
//		for i in 0..<buffer {
//			result.flatValues[i] = lhs.flatValues[i] / rhs
//		}
//		return result
//	}
//	/// Scalar multiplication
//	public static func * (lhs: T, rhs: InlineMatrix<buffer,T>) -> InlineMatrix<buffer,T> {
//		var result = rhs
//		for i in 0..<buffer {
//			result.flatValues[i] = rhs.flatValues[i] * lhs
//		}
//		return result
//	}
//	
//	/// Elementwise multiplication of two matricies
//	public mutating func multiplyElementwise(by other: Self) {
//		for i in 0..<buffer {
//			flatValues[i] = flatValues[i]*other.flatValues[i]
//		}
//	}
//	/// Elementwise multiplication of two matricies
//	public func multiplyingElementwise(by other: Self) -> Self {
//		var newMatrix = self
//		newMatrix.multiplyElementwise(by: other)
//		return newMatrix
//	}
//	
//	/// Perform Matrix Multiplication
//	public static func * <let N: Int, let M: Int>(lhs: InlineMatrix<N,T>, rhs: InlineMatrix<M,T>) throws -> InlineMatrix<buffer,T> {
//		guard lhs.columns == rhs.rows else { throw MatrixError.nonMatchingDimensions }
//		var result = InlineMatrix<buffer,T>(columns: rhs.columns)
//		for i in 0..<lhs.rows {
//			for j in 0..<rhs.columns {
//				for k in 0..<lhs.columns {
//					result.values[i][j] += lhs.values[i][k] * rhs.values[k][j]
//				}
//			}
//		}
//		return result
//	}
//	
//	/// Returns the matrix as a single inline array, only if matrix is a single column or a single row.
//	public func makeInlineVector() throws -> InlineArray<buffer,T> {
//		guard rows == 1 || columns == 1 else { throw MatrixError.nonMatchingDimensions }
//		if rows == 1 {
//			guard columns == buffer else { throw MatrixError.nonMatchingDimensions }
//		}
//		if columns == 1 {
//			guard rows == buffer else { throw MatrixError.nonMatchingDimensions }
//		}
//		return flatValues
//	}
//	/// Returns the matrix as a single array, only if matrix is a single column or a single row.
//	public func makeVectorArray() throws -> Array<T> {
//		guard rows == 1 || columns == 1 else { throw MatrixError.nonMatchingDimensions }
//		if rows == 1 {
//			guard columns == buffer else { throw MatrixError.nonMatchingDimensions }
//		}
//		if columns == 1 {
//			guard rows == buffer else { throw MatrixError.nonMatchingDimensions }
//		}
//		return Array(flatValues)
//	}
//	
//	/// Transposes the matrix in-place (swapping rows for columns)
//	public mutating func transposeSelf() {
//		guard buffer > 0 else { return }
//		for r in 0..<rows {
//			for c in 0..<columns {
//				let i1 = r * columns + c
//				let i2 = c * rows + r
//				flatValues.swapAt(i1, i2)
//			}
//		}
////		flatValues.withUnsafeMutableBufferPointer { pointer in
////			vDSP_mtransD(flatValues, 1, pointer.baseAddress!, 1, vDSP_Length(columns), vDSP_Length(rows))
////		}
//	}
//	/// Returns the transpose of the matrix (swapping rows for columns)
//	public func transpose() -> Self {
//		guard rows > 0 && columns > 0 else { return self }
//		var result = self
//		result.transposeSelf()
//		return result
//	}
//	
//	/// Returns `true` if matrix is symmetric
//	public var isSymmetric: Bool {
//		if isSquare {
//			for r in 0..<rows {
//				for c in 0..<columns {
//					let i1 = r * columns + c
//					let i2 = c * rows + r
//					if flatValues[i1] != flatValues[i2] { return false }
//				}
//			}
//			return true
//		} else { return false }
//	}
//	
//	/// Check if the matrix is square.
//	private var isSquare: Bool { rows == columns }
//	
//	/// Determinant of the matrix (only defined for square matrices).
//	public func determinant() throws -> T {
//		guard isSquare else { throw MatrixError.notSquareMatrix }
//		return InlineMatrix.calculateDeterminant(for: values)
//	}
//	
//	public func cofactorMatrix() throws -> InlineMatrix<buffer, T> {
//		guard isSquare else {
//			throw MatrixError.notSquareMatrix
//		}
//		let n = rows
//		guard n > 1 else { throw MatrixError.indexOutOfRange }
//		if n == 2 {
//			return InlineMatrix([
//				[self[1,1], -self[1,0]],
//				[-self[0,1], self[0,0]]
//			])
//		}
//		var cofactorMatrix = InlineMatrix<buffer, T>(columns: columns)
//		for i in 0..<n {
//			for j in 0..<n {
//				let minor = minorMatrix(removingRow: i, removingColumn: j)
//				let sign: T = ((i + j) % 2 == 0) ? 1 : -1
//				cofactorMatrix[i][j] = sign * InlineMatrix.calculateDeterminant(for: minor)
//			}
//		}
//		return InlineMatrix(cofactorMatrix)
//	}
//	public func adjoint() throws -> InlineMatrix<buffer, T> {
//		// Transpose the cofactor matrix to get the adjoint
//		return try self.cofactorMatrix().transpose()
//	}
//	private func minorMatrix(removingRow rowToRemove: Int, removingColumn columnToRemove: Int) -> [[T]] {
//		values.enumerated().compactMap { i, row in
//			guard i != rowToRemove else { return nil }
//			let filteredRow = row.enumerated().compactMap { j, value in
//				j != columnToRemove ? value : nil
//			}
//			return filteredRow
//		}
//	}
//	
//	public func removing(rows rowsToRemove: [Int], columns columnsToRemove: [Int]) throws -> Matrix<T> {
//		for row in rowsToRemove {
//			if row < 0 || row >= rows { throw MatrixError.indexOutOfRange }
//		}
//		for column in columnsToRemove {
//			if column < 0 || column >= columns { throw MatrixError.indexOutOfRange }
//		}
//		let rowSet = Set(rowsToRemove)
//		let columnSet = Set(columnsToRemove)
//		let filteredValues: [[T]] = values.enumerated().compactMap { (i, row) in
//			guard !rowSet.contains(i) else { return nil }
//			let filteredRow = row.enumerated().compactMap { (j, value) in
//				columnSet.contains(j) ? nil : value
//			}
//			return filteredRow
//		}
//		return Matrix(filteredValues)
//	}
//	
//	/// Recursive helper function to calculate the determinant.
//	private static func calculateDeterminant(for matrix: [[T]]) -> T {
//		let n = matrix.count
//	
//		if n == 1 { // 1x1 matrix
//			return matrix[0][0]
//		}
//		if n == 2 { // 2x2 matrix
//			return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]
//		}
//		if n == 3 { // 3x3 matrix
//			let a = matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1])
//			let b = matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0])
//			let c = matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0])
//			return a - b + c
//		}
//		if n < 6 {
//			// Calculate determinant using cofactor expansion - O(n!)
//			var det: T = 0
//			for i in 0..<n {
//				// Create the submatrix by excluding the first row and the ith column
//				var subMatrix = [[T]]()
//				for j in 1..<matrix.count {
//					var row = [T]()
//					for k in 0..<matrix.count {
//						if k != i {
//							row.append(matrix[j][k])
//						}
//					}
//					subMatrix.append(row)
//				}
//				// Recursive determinant calculation with cofactor expansion
//				let sign: T = (i % 2 == 0) ? 1 : -1
//				det += sign * matrix[0][i] * calculateDeterminant(for: subMatrix)
//			}
//			return det
//		} else {
//			// Calculate determinant using LU decomposition - O(n^3)
//			do {
//				let (L, U, _, rowSwaps) = try Matrix.luDecomposition(of: matrix)
//				return calculateDeterminant(fromL: L, withU: U, rowSwaps: rowSwaps)
//			} catch {
//				if let matrixError = error as? MatrixError {
//					if matrixError == MatrixError.singularMatrix { return 0 }
//				}
//				fatalError(error.localizedDescription)
//			}
//		}
//	}
//	/// Calculate the determinant from LU Decomposition
//	private static func calculateDeterminant(fromL L: [[T]], withU U: [[T]], rowSwaps: Int) -> T {
//		//let detL: T = 1
//		var detU: T = 1
//		// Product of diagonal elements
//		for i in 0..<L.count {
//			//detL *= L[i][i] // Diagonal elements are all 1.
//			detU *= U[i][i]
//		}
//		//let det = detL*detU
//		// Adjust sign for number of row swaps
//		return (rowSwaps % 2 == 0) ? detU : -detU
//	}
//	
//	public static func solve(A: InlineMatrix<buffer, T>, b: InlineMatrix<buffer, T>) throws -> InlineMatrix<buffer, T> {
//		guard A.isSquare else { throw MatrixError.notSquareMatrix }
//		guard b.columns == 1 && b.rows == A.rows else { throw MatrixError.nonMatchingDimensions }
//		let n = A.rows
//		if n < 6 {
//			let inverseA = try A.inverse()
//			return try inverseA * b
//		}
//		
//		let (L, U, P, _) = try InlineMatrix.luDecomposition(of: A.values)
//		
//		let x = try Matrix.solve(L: L, U: U, P: P, b: Array(b.flatValues))
//		return InlineMatrix(x, isRow: false)
//	}
//	public static func solve(L: [[T]], U: [[T]], P: [Int], b: [T]) throws -> [T] {
//		guard L.count == L.first?.count && U.count == U.first?.count else {
//			throw MatrixError.notSquareMatrix
//		}
//		guard b.count == L.count && L.count == U.count && P.count == L.count else {
//			throw MatrixError.nonMatchingDimensions
//		}
//		let n = L.count
//		
//		// Helper to apply permutation to identity column
//		func permute(_ e: [T]) -> [T] {
//			var result = Array<T>(repeating: 0.0, count: n)
//			for i in 0..<n {
//				result[i] = e[P[i]]
//			}
//			return result
//		}
//		
//		let Pb = permute(b)
//		
//		// Forward substitution: solve L * y = Pb
//		var y = Array<T>(repeating: 0.0, count: n)
//		for j in 0..<n {
//			y[j] = Pb[j]
//			for k in 0..<j {
//				y[j] -= L[j][k] * y[k]
//			}
//			// L diagonal is 1, so divide unnecessary
//		}
//		
//		// Backward substitution: solve U * x = y
//		var x = Array<T>(repeating: 0.0, count: n)
//		for j in stride(from: n - 1, through: 0, by: -1) {
//			x[j] = y[j]
//			for k in (j + 1)..<n {
//				x[j] -= U[j][k] * x[k]
//			}
//			guard U[j][j] != 0 else { throw MatrixError.singularMatrix } // Singular matrix
//			x[j] /= U[j][j]
//		}
//		return x
//	}
//	
//	/// Determinant of the matrix (only defined for square matrices).
//	public func inverse() throws -> InlineMatrix<buffer, T> {
//		guard isSquare else { throw MatrixError.notSquareMatrix }
//		let n = rows
//		if n < 6 {
//			let det = InlineMatrix.calculateDeterminant(for: values)
//			if det == 0 { throw MatrixError.singularMatrix }
//			
//			if n == 1 { return InlineMatrix([[(1/self[0,0])]]) }
//			if n == 3 {
//				let m = values
//
//				let a = m[0][0], b = m[0][1], c = m[0][2]
//				let d = m[1][0], e = m[1][1], f = m[1][2]
//				let g = m[2][0], h = m[2][1], i = m[2][2]
//
//				let A = e * i - f * h
//				let B = -(d * i - f * g)
//				let C = d * h - e * g
//				let D = -(b * i - c * h)
//				let E = a * i - c * g
//				let F = -(a * h - b * g)
//				let G = b * f - c * e
//				let H = -(a * f - c * d)
//				let I = a * e - b * d
//				
//				//let det = a * A + b * B + c * C
//				let resultValues: [[T]] = [
//					[A, D, G],
//					[B, E, H],
//					[C, F, I]
//				]
//				return InlineMatrix(resultValues)/det
//			}
//			// Use Cofactor Expansion for small matrices (n! < n^3).
//			return (1/det)*(try! self.adjoint())
//		}
//		
//		// TODO: Implement Cholesky decomposition for symmetric positive definite matrices
//		// Use LU Decomposition
//		let (L, U, P, _) = try InlineMatrix.luDecomposition(of: self.values)
//		var invCols: [[T]] = []
//		
//		// Helper to apply permutation to identity column
//		func permute(_ e: [T]) -> [T] {
//			var result = Array<T>(repeating: 0.0, count: n)
//			for i in 0..<n {
//				result[i] = e[P[i]]
//			}
//			return result
//		}
//		
//		// Solve for each column of the inverse
//		for i in 0..<n {
//			// Unit vector (column of identity)
//			var e = Array<T>(repeating: 0.0, count: n)
//			e[i] = 1.0
//			let b = permute(e)
//			
//			// Forward substitution: solve L * y = b
//			var y = Array<T>(repeating: 0.0, count: n)
//			for j in 0..<n {
//				y[j] = b[j]
//				for k in 0..<j {
//					y[j] -= L[j][k] * y[k]
//				}
//				// L diagonal is 1, so divide unnecessary
//			}
//			
//			// Backward substitution: solve U * x = y
//			var x = Array<T>(repeating: 0.0, count: n)
//			for j in stride(from: n - 1, through: 0, by: -1) {
//				x[j] = y[j]
//				for k in (j + 1)..<n {
//					x[j] -= U[j][k] * x[k]
//				}
//				guard U[j][j] != 0 else { throw MatrixError.singularMatrix } // Singular matrix
//				x[j] /= U[j][j]
//			}
//			
//			invCols.append(x)
//		}
//		
//		// Transpose columns to rows
//		let inverseValues = (0..<n).map { i in
//			invCols.map { $0[i] }
//		}
//		return InlineMatrix(inverseValues)
//
//		// Perform Gauss-Jordan elimination
////		var augmented = values
////
////		// Append identity matrix to the right of the original
////		for i in 0..<size {
////			augmented[i] += (0..<size).map { $0 == i ? 1 : 0 }.map(T.init)
////		}
////
////		for i in 0..<size {
////			// Find pivot
////			var pivotRow = i
////			for j in i+1..<size where abs(augmented[j][i]) > abs(augmented[pivotRow][i]) {
////				pivotRow = j
////			}
////
////			// If pivot is zero, matrix is singular
////			if augmented[pivotRow][i] == 0 {
////				throw MatrixError.singularMatrix
////			}
////
////			// Swap rows if needed
////			if i != pivotRow {
////				augmented.swapAt(i, pivotRow)
////			}
////
////			// Normalize pivot row
////			let pivot = augmented[i][i]
////			for j in 0..<2 * size {
////				augmented[i][j] /= pivot
////			}
////
////			// Eliminate other rows
////			for k in 0..<size where k != i {
////				let factor = augmented[k][i]
////				for j in 0..<2 * size {
////					augmented[k][j] -= factor * augmented[i][j]
////				}
////			}
////		}
////
////		// Extract right half as the inverse
////		let inverseValues = augmented.map { Array($0[size..<(2*size)]) }
////		return Matrix(inverseValues)
//	}
//	
//	/// Returns (L, U, P, swapCount) for LU decomposition with partial pivoting.
//	/// - Returns: Optional tuple (L, U, P, swapCount), or nil if decomposition fails.
//	public static func luDecomposition(of matrix: [[T]], tolerance: T = 0.000001) throws -> (L: [[T]], U: [[T]], P: [Int], swapCount: Int) {
//		guard matrix.count == matrix.first?.count else { throw MatrixError.notSquareMatrix }
//		let n = matrix.count
//		var A = matrix
//		var L = Array(repeating: Array<T>(repeating: 0.0, count: n), count: n)
//		var U = Array(repeating: Array<T>(repeating: 0.0, count: n), count: n)
//		var P: Array<Int> = Array(0..<n)
//		var swapCount = 0
//		
//		let swapsAreRequired: Bool = {
//			for i in 0..<n {
//				if abs(A[i][i]) < tolerance { return true }
//			}
//			return false
//		}()
//		
//		if swapsAreRequired {
//			for i in 0..<n {
//				// Partial Pivoting
//				var maxRow = i
//				var maxVal: T = 0.0
//				for k in i..<n {
//					let kABS = abs(A[k][i])
//					if kABS > maxVal {
//						maxVal = kABS
//						maxRow = k
//					}
//				}
//				
//				// Zero pivot => singular (or near singular) matrix
//				if maxVal < tolerance {
//					throw MatrixError.singularMatrix
//				}
//				
//				// Swap rows in A and P
//				if maxRow != i {
//					A.swapAt(i, maxRow) // Pivot A
//					P.swapAt(i, maxRow) // Pivot P
//					swapCount += 1
//				}
//			}
//		}
//
//		// Doolittle algorithm + Partial Pivoting
//		for i in 0..<n {
//			
//			// Upper Triangular
//			for j in i..<n {
//				// Summation of L(i, k) * U(k, j)
//				var sum: T = 0.0
//				for k in 0..<i {
//					sum += L[i][k] * U[k][j]
//				}
//				// Evaluating U(i, j)
//				U[i][j] = A[i][j] - sum
//			}
//
//			
//			// Lower Triangular
//			L[i][i] = 1.0 // unit diagonal
//			for j in (i+1)..<n {
//				// Summation of L(j, k) * U(k, i)
//				var sum: T = 0.0
//				for k in 0..<i {
//					sum += L[j][k] * U[k][i]
//				}
//				// Evaluating L(j, i)
//				if U[i][i] == 0.0 { throw MatrixError.singularMatrix } // Singular matrix
//				L[j][i] = (A[j][i] - sum) / U[i][i]
//			}
//		}
//		return (L, U, P, swapCount)
//	}
//	/// Builds a permutation matrix from a row permutation array.
//	/// Each index `i` in `permutation` maps original row `i` to `permutation[i]`.
//	public static func permutationMatrix(from permutation: [Int]) -> InlineMatrix<buffer, T> {
//		let n = permutation.count
//		var values = Array(repeating: Array<T>(repeating: 0.0, count: n), count: n)
//		for i in 0..<n {
//			values[permutation[i]][i] = 1.0
//		}
//		return InlineMatrix(values)
//	}
//	
//	// UNTESTED - written by AI - LU decomposition is about 2x faster than QR usually for solving matrices, but QR is best for eigenvalues and eigenvectors.
//	private func qrDecompose(iterations: Int = 10) -> (Q: [[T]], R: [[T]]) {
//		var A = self.values
//		let n = A.count
//		let m = A[0].count
//		var Q: [[T]] = Array(repeating: Array(repeating: 0.0, count: m), count: n)
//		var R: [[T]] = Array(repeating: Array(repeating: 0.0, count: m), count: m)
//
//		for k in 0..<m {
//			var norm: T = 0.0
//			for i in 0..<n {
//				norm += A[i][k] * A[i][k]
//			}
//			norm = sqrt(norm)
//			for i in 0..<n {
//				Q[i][k] = A[i][k] / norm
//			}
//			R[k][k] = norm
//
//			for j in (k+1)..<m {
//				var dot: T = 0.0
//				for i in 0..<n {
//					dot += Q[i][k] * A[i][j]
//				}
//				R[k][j] = dot
//				for i in 0..<n {
//					A[i][j] -= Q[i][k] * dot
//				}
//			}
//		}
//		return (Q, R)
//	}
//
//	/// Display row count x column count
//	public var sizeDescription: String {
//		"\(rows)x\(columns)"
//	}
//	
//	/// Convert matrix to a string representation
//	public var description: String {
//		var result = "\(rows)x\(columns)\n"
//		for i in 0..<rows {
//			result += "["
//			for j in 0..<columns {
//				result += "\(self[i,j])"
//				if j < columns - 1 {
//					result += ", "
//				}
//			}
//			result += "]\n"
//		}
//		return result
//	}
//	
//	/// Check if indices are valid
//	private func isValidIndex(row: Int, col: Int) -> Bool {
//		return row >= 0 && row < rows && col >= 0 && col < columns
//	}
//	/// Fetches or updates a value in the matrix. Will crash if out-of-range.
//	public subscript(r: Int, c: Int) -> T {
//		get {
//			guard isValidIndex(row: r, col: c) else {
//				fatalError("Matrix index out of bounds: (\(r), \(c))")
//			}
//			return flatValues[r * columns + c]
//		}
//		set(newValue) {
//			guard isValidIndex(row: r, col: c) else {
//				fatalError("Matrix index out of bounds: (\(r), \(c))")
//			}
//			flatValues[r * columns + c] = newValue
//		}
//	}
//	
//	// MARK: Surge
//	// Surge MIT License:
//
//	//Copyright © 2014-2019 the Surge contributors
//	//
//	//Permission is hereby granted, free of charge, to any person obtaining a copy
//	//of this software and associated documentation files (the "Software"), to deal
//	//in the Software without restriction, including without limitation the rights
//	//to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
//	//copies of the Software, and to permit persons to whom the Software is
//	//furnished to do so, subject to the following conditions:
//	//
//	//The above copyright notice and this permission notice shall be included in
//	//all copies or substantial portions of the Software.
//	//
//	//THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//	//IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//	//FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
//	//AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//	//LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
//	//OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
//	//THE SOFTWARE.
//	
//	public func eigenvalues(maxIterations: Int = 100, tolerance: Double = 1e-10) throws -> [(x: T, i: T)] {
//		guard isSquare else { throw MatrixError.notSquareMatrix }
//		if rows == 1 {
//			return self[0,0].isNaN ? [] : [((self[0,0]), T(0))]
//		}
//		if rows == 2 {
//			let a = self[0,0]
//			let b = self[0,1]
//			let c = self[1,0]
//			let d = self[1,1]
//
//			let trace = a + d
//			let determinant = a * d - b * c
//			let discriminant = trace * trace - 4 * determinant
//
//			if discriminant >= 0 {
//				let sqrtDiscriminant = discriminant.squareRoot()
//				let lambda1 = (trace + sqrtDiscriminant) / 2
//				let lambda2 = (trace - sqrtDiscriminant) / 2
//				return [(lambda1,0.0), (lambda2,0.0)]
//			}
//		}
//		return try eigenDecompose(computeEigenVectors: false).eigenValues
//	}
//	
//	/// Decomposes a square matrix into its eigenvalues and left and right eigenvectors.
//	/// The decomposition may result in complex numbers, represented by (T, T), which
//	///   are the (real, imaginary) parts of the complex number.
//	/// - Returns: a struct with the eigen values and left and right eigen vectors using (T, T)
//	///   to represent a complex number.
//	func eigenDecompose(computeEigenVectors: Bool) throws -> MatrixEigenDecompositionResult<T> {
//		let input = Matrix<Double>(self)
//		let decomposition = try eigenDecomposeOLD(input, computeEigenVectors: computeEigenVectors)
//		
////		guard rows == columns else {
////			throw MatrixError.notSquareMatrix
////		}
////		let (Q,R) = self.qrDecompose()
////		let lambdas: [(x: T, i: T)] = []
////		var rightEigenVectors: [[(T, T)]] = [[]]
////		if computeEigenVectors  {
////			for lambda in lambdas {
////				var A = try self-Matrix<T>.identity(size: rows)*lambda.x
////				let realVector = try Matrix.solve(A: A, b: Matrix(Array<T>(repeating: 0, count: rows), isRow: false))
////				A = try self-Matrix<T>.identity(size: rows)*lambda.i
////				let imaginaryVector = try Matrix.solve(A: A, b: Matrix(Array<T>(repeating: 0, count: rows), isRow: false))
////				var eigenVector: [(T, T)] = []
////				for i in 0..<rows {
////					eigenVector.append((realVector.flatValues[i],imaginaryVector.flatValues[i]))
////				}
////				rightEigenVectors.append(eigenVector)
////			}
////		}
//		
//		return MatrixEigenDecompositionResult<T>(
//			eigenValues: decomposition.eigenValues.map { (T($0.0), T($0.1)) },
//			leftEigenVectors: decomposition.leftEigenVectors.map { $0.map { (T($0.0), T($0.1)) } },
//			rightEigenVectors: decomposition.rightEigenVectors.map { $0.map { (T($0.0), T($0.1)) } }
//		)
//	}
//
//	/// Decomposes a square matrix into its eigenvalues and left and right eigenvectors.
//	/// The decomposition may result in complex numbers, represented by (Double, Double), which
//	///   are the (real, imaginary) parts of the complex number.
//	/// - Parameters:
//	///   - lhs: a square matrix
//	/// - Returns: a struct with the eigen values and left and right eigen vectors using (Double, Double)
//	///   to represent a complex number.
//	func eigenDecomposeOLD(_ lhs: Matrix<Double>, computeEigenVectors: Bool) throws -> MatrixEigenDecompositionResult<Double> {
//		
//		guard lhs.rows == lhs.columns else {
//			throw MatrixError.notSquareMatrix
//		}
//
//		// dgeev_ needs column-major matrices, so transpose 'lhs'.
//		var matrixGrid: [__CLPK_doublereal] = lhs.transpose().flatValues
//		var matrixRowCount = __CLPK_integer(lhs.rows)
//		let matrixColCount = matrixRowCount
//		var eigenValueCount = matrixRowCount
//		var leftEigenVectorCount = matrixRowCount
//		var rightEigenVectorCount = matrixRowCount
//
//		var workspaceQuery: Double = 0.0
//		var workspaceSize = __CLPK_integer(-1)
//		var error: __CLPK_integer = 0
//
//		var eigenValueRealParts = [Double](repeating: 0, count: Int(eigenValueCount))
//		var eigenValueImaginaryParts = [Double](repeating: 0, count: Int(eigenValueCount))
//		var leftEigenVectorWork = [Double](repeating: 0, count: Int(leftEigenVectorCount * matrixColCount))
//		var rightEigenVectorWork = [Double](repeating: 0, count: Int(rightEigenVectorCount * matrixColCount))
//
//		var decompositionJobVL: [CChar]
//		var decompositionJobVR: [CChar]
//		
//		if computeEigenVectors  {
//			decompositionJobVL = [0x56, 0x00] // "V" (compute)
//			decompositionJobVR = [0x56, 0x00] // "V" (compute)
//		} else {
//			decompositionJobVL = Array("N".utf8CString) // "N" (do not compute)
//			decompositionJobVR = Array("N".utf8CString) // "N" (do not compute)
//		}
//
//		// Call dgeev to find out how much workspace to allocate
//		dgeev_(&decompositionJobVL, &decompositionJobVR, &matrixRowCount, &matrixGrid, &eigenValueCount, &eigenValueRealParts, &eigenValueImaginaryParts, &leftEigenVectorWork, &leftEigenVectorCount, &rightEigenVectorWork, &rightEigenVectorCount, &workspaceQuery, &workspaceSize, &error)
//		if error != 0 {
//			throw MatrixError.matrixNotDecomposable
//		}
//
//		// Allocate the workspace and call dgeev again to do the actual decomposition
//		var workspace = [Double](repeating: 0.0, count: Int(workspaceQuery))
//		workspaceSize = __CLPK_integer(workspaceQuery)
//		dgeev_(&decompositionJobVL, &decompositionJobVR, &matrixRowCount, &matrixGrid, &eigenValueCount, &eigenValueRealParts, &eigenValueImaginaryParts, &leftEigenVectorWork, &leftEigenVectorCount, &rightEigenVectorWork, &rightEigenVectorCount, &workspace, &workspaceSize, &error)
//		if error != 0 {
//			throw MatrixError.matrixNotDecomposable
//		}
//
//		return MatrixEigenDecompositionResult<Double>(rowCount: lhs.rows, eigenValueRealParts: eigenValueRealParts, eigenValueImaginaryParts: eigenValueImaginaryParts, leftEigenVectorWork: leftEigenVectorWork, rightEigenVectorWork: rightEigenVectorWork)
//	}
//	
//	/// Holds the result of eigendecomposition. The (Scalar, Scalar) used
//	/// in the property types represents a complex number with (real, imaginary) parts.
//	struct MatrixEigenDecompositionResult<Scalar: BinaryFloatingPoint> {
//		public let eigenValues: [(Scalar, Scalar)]
//		public let leftEigenVectors: [[(Scalar, Scalar)]]
//		public let rightEigenVectors: [[(Scalar, Scalar)]]
//		
//		public init(eigenValues: [(Scalar, Scalar)], leftEigenVectors: [[(Scalar, Scalar)]], rightEigenVectors: [[(Scalar, Scalar)]]) {
//			self.eigenValues = eigenValues
//			self.leftEigenVectors = leftEigenVectors
//			self.rightEigenVectors = rightEigenVectors
//		}
//		
//		public init(rowCount: Int, eigenValueRealParts: [Scalar], eigenValueImaginaryParts: [Scalar], leftEigenVectorWork: [Scalar], rightEigenVectorWork: [Scalar]) {
//			// The eigenvalues are an array of (real, imaginary) results from dgeev
//			self.eigenValues = Array(zip(eigenValueRealParts, eigenValueImaginaryParts))
//			
//			// Build the left and right eigenvectors
//			let emptyVector = [(Scalar, Scalar)](repeating: (0.0, 0.0), count: rowCount)
//			var leftEigenVectors = [[(Scalar, Scalar)]](repeating: emptyVector, count: rowCount)
//			MatrixEigenDecompositionResult.buildEigenVector(eigenValueImaginaryParts: eigenValueImaginaryParts, eigenVectorWork: leftEigenVectorWork, result: &leftEigenVectors)
//			
//			var rightEigenVectors = [[(Scalar, Scalar)]](repeating: emptyVector, count: rowCount)
//			MatrixEigenDecompositionResult.buildEigenVector(eigenValueImaginaryParts: eigenValueImaginaryParts, eigenVectorWork: rightEigenVectorWork, result: &rightEigenVectors)
//			
//			self.leftEigenVectors = leftEigenVectors
//			self.rightEigenVectors = rightEigenVectors
//		}
//		
//		// Convert the result of dgeev into an array of complex numbers
//		// See Intel's documentation on column-major results for sample code that this is based on:
//		// https://software.intel.com/sites/products/documentation/doclib/mkl_sa/11/mkl_lapack_examples/dgeev.htm
//		private static func buildEigenVector<S>(eigenValueImaginaryParts: [S], eigenVectorWork: [S], result: inout [[(S, S)]]) where S: FloatingPoint & ExpressibleByFloatLiteral {
//			// row and col count are the same because result must be square.
//			let rowColCount = result.count
//
//			for row in 0..<rowColCount {
//				var col = 0
//				while col < rowColCount {
//					if eigenValueImaginaryParts[col] == 0.0 {
//						// v is column-major
//						result[row][col] = (eigenVectorWork[row + rowColCount * col], 0.0)
//						col += 1
//					} else {
//						// v is column-major
//						result[row][col] = (eigenVectorWork[row + col * rowColCount], eigenVectorWork[row + rowColCount * (col + 1)])
//						result[row][col + 1] = (eigenVectorWork[row + col * rowColCount], -eigenVectorWork[row + rowColCount * (col + 1)])
//						col += 2
//					}
//				}
//			}
//		}
//
//	}
//}

//@available(iOS 26.0.0, macOS 26.0.0, *)
//nonisolated private struct InlineMatrix<let rows: Int, let columns: Int, let bufferSize: Int, T: BinaryFloatingPoint>: Equatable, CustomStringConvertible, Copyable {
//	
////	public var rows: Int { rows }
////	public var columns: Int { columns }
//	public var values: [[T]] {
//		get {
//			var result = [[T]]()
//			for i in 0..<bufferSize {
//				let r = i / columns
//				let c = i % columns
//				result[r][c] = flatValues[i]
//			}
//			return result
//		}
//		set {
//			for i in 0..<bufferSize {
//				let r = i / columns
//				let c = i % columns
//				flatValues[i] = newValue[r][c]
//			}
//		}
//	}
//	public private(set)var flatValues: InlineArray<bufferSize, T>
//	
//	static func == (lhs: InlineMatrix<rows, columns, bufferSize, T>, rhs: InlineMatrix<rows, columns, bufferSize, T>) -> Bool {
//		lhs.sizeDescription == rhs.sizeDescription && lhs.flatValues == rhs.flatValues
//	}
//	func isIdentical(to rhs: InlineMatrix<rows, columns, bufferSize, T>) -> Bool {
//		self.flatValues.span.isIdentical(to: rhs.flatValues.span)
//	}
//	
//	public init<V>(_ matrix: InlineMatrix<rows, columns, bufferSize, V>) where V: BinaryFloatingPoint {
//		self.flatValues = InlineArray<bufferSize, T>() { T(matrix.flatValues[$0]) }
//	}
//	/// Creates a matrix where the outer array is for rows and the inner array is for columns.
//	public init(_ values: [[T]]) {
//		let tempFlatValues = values.flatMap(\.self)
//		self.flatValues = InlineArray() {
//			tempFlatValues[$0]
//		}
//	}
//	public init(_ values: [T]) {
//		self.flatValues = InlineArray() { values[$0] }
//	}
//	public init(_ values: InlineArray<bufferSize, T>) {
//		self.flatValues = values
//	}
//	
//	/// Initialize a matrix of size n x m with a default value
//	public init(repeating value: T = 0) {
//		self.flatValues = InlineArray<bufferSize, T>(repeating: value)
//	}
//	
//	public init(_ closure: (_ row: Int, _ column: Int) throws -> T) rethrows {
//		self.flatValues = try InlineArray<bufferSize, T>() { i in
//			let r = i / columns
//			let c = i % columns
//			return try closure(r, c)
//		}
//	}
//	public init(_ closure: (_ i: Int) throws -> T) rethrows {
//		self.flatValues = try InlineArray<bufferSize, T>() { i in
//			return try closure(i)
//		}
//	}
//
//	public static func identity() -> InlineMatrix<rows, columns, bufferSize, T> {
//		return InlineMatrix<rows, columns, bufferSize, T>() { i in
//			if i / columns == i % columns { 1 } else { 0 }
//		}
//	}
//	public static func diagonal(repeating value: T) -> InlineMatrix<rows, columns, bufferSize, T> {
//		let n = min(rows,columns)
//		return InlineMatrix<rows, columns, bufferSize, T>() { i in
//			if i / n == i % n { value } else { 0 }
//		}
//	}
//	public static func diagonal<C>(scalars: C) -> InlineMatrix<rows, columns, bufferSize, T> where C: Collection, C.Element == T {
//		let n = min(rows,columns)
//		var iScalar: C.Index = scalars.index(scalars.startIndex, offsetBy: -1)
//		return InlineMatrix<rows, columns, bufferSize, T>() { i in
//			if i / n == i % n {
//				iScalar = scalars.index(after: iScalar)
//				return scalars[iScalar]
//			} else { return 0 }
//		}
//	}
//
//	public init(_ simdMatrix: simd_double2x2) where T == Double, rows == 2, columns == 2 {
//		let m = simdMatrix.transpose
//		let temp = [
//			m.columns.0.x, m.columns.0.y,
//			m.columns.1.x, m.columns.1.y,
//		]
//		self.flatValues = InlineArray() { temp[$0] }
//	}
//	public init(_ simdMatrix: simd_double3x3) where T == Double, rows == 3, columns == 3 {
//		let m = simdMatrix.transpose
//		let temp = [
//			m.columns.0.x, m.columns.0.y, m.columns.0.z,
//			m.columns.1.x, m.columns.1.y, m.columns.1.z,
//			m.columns.2.x, m.columns.2.y, m.columns.2.z,
//		]
//		self.flatValues = InlineArray() { temp[$0] }
//	}
//	public init(_ simdMatrix: simd_double4x4) where T == Double, rows == 4, columns == 4 {
//		let m = simdMatrix.transpose
//		let temp = [
//			m.columns.0.x, m.columns.0.y, m.columns.0.z, m.columns.0.w,
//			m.columns.1.x, m.columns.1.y, m.columns.1.z, m.columns.1.w,
//			m.columns.2.x, m.columns.2.y, m.columns.2.z, m.columns.2.w,
//			m.columns.3.x, m.columns.3.y, m.columns.3.z, m.columns.3.w,
//		]
//		self.flatValues = InlineArray() { temp[$0] }
//	}
//	
//	/// Check if rows and columns of each matrix match.
//	public static func areDimensionsEqual<let N1: Int, let N2: Int, let M1: Int, let M2: Int, let S1: Int, let S2: Int>(lhs: InlineMatrix<N1,M1,S1,T>, rhs: InlineMatrix<N2,M2,S2,T>) -> Bool {
//		N1 == N2 && M1 == M2 && S1 == S2
//	}
//	
//	public subscript(row row: Int) -> InlineArray<columns, T> {
//		get {
//			assert(row < rows)
//			let startIndex = row * columns
//			return InlineArray<columns, T>() { c in
//				flatValues[startIndex+c]
//			}
//		}
//		set {
//			assert(row < rows)
//			assert(newValue.count == columns)
//			let startIndex = row * columns
//			for c in 0..<newValue.count {
//				flatValues[c+startIndex] = newValue[c]
//			}
//		}
//	}
//	public subscript(column column: Int) -> InlineArray<rows, T> {
//		get {
//			return InlineArray<rows, T>() { r in
//				let index = r * columns + column
//				return self.flatValues[index]
//			}
//		}
//		set {
//			assert(column < columns)
//			assert(newValue.count == rows)
//			for r in 0..<rows {
//				let index = r * columns + column
//				flatValues[index] = newValue[r]
//			}
//		}
//	}
//	
//	/// Returns `true` if both matrices are approximately equivalent to 10^-4
//	public static func ~= (lhs: InlineMatrix<rows, columns, bufferSize, T>, rhs: InlineMatrix<rows, columns, bufferSize,T>) -> Bool {
//		guard InlineMatrix.areDimensionsEqual(lhs: lhs, rhs: rhs) else { return false }
//		let result = try! lhs - rhs
//		for i in 0..<bufferSize {
//			if abs(result.flatValues[i]) > 0.000001 {
//				return false
//			}
//		}
//		return true
//	}
//	
//	/// Elementwise addition
//	public static func + (lhs: InlineMatrix<rows, columns, bufferSize,T>, rhs: InlineMatrix<rows, columns, bufferSize,T>) throws -> InlineMatrix<rows, columns, bufferSize,T> {
//		guard InlineMatrix.areDimensionsEqual(lhs: lhs, rhs: rhs) else { throw MatrixError.nonMatchingDimensions }
//		var result = lhs
//		for i in 0..<bufferSize {
//			result.flatValues[i] = lhs.flatValues[i] - rhs.flatValues[i]
//		}
//		return result
//	}
//	/// Elementwise subtraction
//	public static func - (lhs: InlineMatrix<rows, columns, bufferSize,T>, rhs: InlineMatrix<rows, columns, bufferSize,T>) throws -> InlineMatrix<rows, columns, bufferSize,T> {
//		guard InlineMatrix.areDimensionsEqual(lhs: lhs, rhs: rhs) else { throw MatrixError.nonMatchingDimensions }
//		var result = lhs
//		for i in 0..<bufferSize {
//			result.flatValues[i] = lhs.flatValues[i] - rhs.flatValues[i]
//		}
//		return result
//	}
//	/// Scalar multiplication
//	public static func * (lhs: InlineMatrix<rows, columns, bufferSize,T>, rhs: T) -> InlineMatrix<rows, columns, bufferSize,T> {
//		var result = lhs
//		for i in 0..<bufferSize {
//			result.flatValues[i] = lhs.flatValues[i] * rhs
//		}
//		return result
//	}
//	/// Scalar division
//	public static func / (lhs: InlineMatrix<rows, columns, bufferSize,T>, rhs: T) -> InlineMatrix<rows, columns, bufferSize,T> {
//		var result = lhs
//		for i in 0..<bufferSize {
//			result.flatValues[i] = lhs.flatValues[i] / rhs
//		}
//		return result
//	}
//	/// Scalar multiplication
//	public static func * (lhs: T, rhs: InlineMatrix<rows, columns, bufferSize,T>) -> InlineMatrix<rows, columns, bufferSize,T> {
//		var result = rhs
//		for i in 0..<bufferSize {
//			result.flatValues[i] = rhs.flatValues[i] * lhs
//		}
//		return result
//	}
//	/// Elementwise multiplication of two matricies
//	public mutating func multiplyElementwise(by other: Self) {
//		for i in 0..<bufferSize {
//			flatValues[i] = flatValues[i]*other.flatValues[i]
//		}
//	}
//	/// Elementwise multiplication of two matricies
//	public func multiplyingElementwise(by other: Self) -> Self {
//		var newMatrix = self
//		newMatrix.multiplyElementwise(by: other)
//		return newMatrix
//	}
//	/// Perform Matrix Multiplication
//	public static func * <let P: Int, let N: Int, let M: Int, let S1: Int, let S2: Int, let S3: Int>(lhs: InlineMatrix<P,N,S1,T>, rhs: InlineMatrix<N,M,S2,T>) -> InlineMatrix<P, M, S3,T> {
//		var result = InlineMatrix<P, M, S3, T>()
//		for i in 0..<P { // lhs rows
//			for j in 0..<M { // rhs columns
//				for k in 0..<N { // lhs colmuns
//					result[i,j] += lhs[i,k] * rhs[k,j]
//				}
//			}
//		}
//		return result
//	}
//	/// Perform Matrix Multiplication
//	public func multiplying<let M: Int, let S: Int, let S2: Int>(by rhs: InlineMatrix<columns,M,S,T>) -> InlineMatrix<rows, M, S2, T> {
//		var result = InlineMatrix<rows, M, S2, T>()
//		for i in 0..<rows { // lhs rows
//			for j in 0..<M { // rhs columns
//				for k in 0..<columns { // lhs colmuns
//					result[i,j] += self[i,k] * rhs[k,j]
//				}
//			}
//		}
//		return result
//	}
//	
//	/// Returns the matrix as a single inline array, only if matrix is a single column or a single row.
//	public func makeInlineVector() throws -> InlineArray<bufferSize,T> {
//		guard rows == 1 || columns == 1 else { throw MatrixError.nonMatchingDimensions }
//		if rows == 1 {
//			guard columns == bufferSize else { throw MatrixError.nonMatchingDimensions }
//		}
//		if columns == 1 {
//			guard rows == bufferSize else { throw MatrixError.nonMatchingDimensions }
//		}
//		return flatValues
//	}
//	/// Returns the matrix as a single array, only if matrix is a single column or a single row.
//	public func makeVectorArray() throws -> Array<T> {
//		return try Array(makeInlineVector())
//	}
//	
//	/// Returns the transpose of the matrix (swapping rows for columns)
//	public func transpose() -> InlineMatrix<columns, rows, bufferSize, T> {
//		return InlineMatrix<columns, rows, bufferSize, T> { row, column in
//			self[column, row]
//		}
//	}
//	
//	/// Returns `true` if matrix is symmetric
//	public var isSymmetric: Bool {
//		if isSquare {
//			for r in 0..<rows {
//				for c in 0..<columns {
//					let i1 = r * columns + c
//					let i2 = c * rows + r
//					if flatValues[i1] != flatValues[i2] { return false }
//				}
//			}
//			return true
//		} else { return false }
//	}
//	
//	/// Check if the matrix is square.
//	private var isSquare: Bool { rows == columns }
//	
//	/// Determinant of the matrix (only defined for square matrices).
//	public func determinant() throws -> T where rows == columns {
//		if rows == 1 { // 1x1 matrix
//			return self[0,0]
//		}
//		if rows == 2 { // 2x2 matrix
//			return self[0,0] * self[1,1] - self[0,1] * self[1,0]
//		}
//		if rows == 3 { // 3x3 matrix
//			let a = self[0,0] * (self[1,1] * self[2,2] - self[1,2] * self[2,1])
//			let b = self[0,1] * (self[1,0] * self[2,2] - self[1,2] * self[2,0])
//			let c = self[0,2] * (self[1,0] * self[2,1] - self[1,1] * self[2,0])
//			return a - b + c
//		}
//		if rows < 6 {
//			// Calculate determinant using cofactor expansion - O(n!)
//			var det: T = 0
//			for i in 0..<rows {
//				// Create the submatrix by excluding the first row and the ith column
//				var subMatrix = [[T]]()
//				for j in 1..<rows {
//					var row = [T]()
//					for k in 0..<rows {
//						if k != i {
//							row.append(self[j,k])
//						}
//					}
//					subMatrix.append(row)
//				}
//				// Recursive determinant calculation with cofactor expansion
//				let sign: T = (i % 2 == 0) ? 1 : -1
//				det += sign * self[0,i] * InlineMatrix.calculateDeterminant(for: subMatrix)
//			}
//			return det
//		} else {
//			// Calculate determinant using LU decomposition - O(n^3)
//			do {
//				let (L, U, _, rowSwaps) = try Matrix.luDecomposition(of: values)
//				return InlineMatrix.calculateDeterminant(fromL: L, withU: U, rowSwaps: rowSwaps)
//			} catch {
//				if let matrixError = error as? MatrixError {
//					if matrixError == MatrixError.singularMatrix { return 0 }
//				}
//				fatalError(error.localizedDescription)
//			}
//		}
//	}
//	
//	public func cofactorMatrix() throws -> InlineMatrix<rows, columns, bufferSize, T> where rows == columns {
//		guard rows > 1 else { throw MatrixError.indexOutOfRange }
//		if rows == 2 {
//			return InlineMatrix([
//				self[1,1], -self[1,0],
//				-self[0,1], self[0,0]
//			])
//		}
//		var cofactorMatrix = InlineMatrix<rows, columns, bufferSize, T>() { r, c in
//			let minor = minorMatrix(removingRow: r, removingColumn: c)
//			let sign: T = ((r + c) % 2 == 0) ? 1 : -1
//			cofactorMatrix[r,c] = sign * InlineMatrix.calculateDeterminant(for: minor)
//		}
//		return InlineMatrix(cofactorMatrix)
//	}
//	public func adjoint() throws -> InlineMatrix<columns, rows, bufferSize, T> where rows == columns {
//		// Transpose the cofactor matrix to get the adjoint
//		return try self.cofactorMatrix().transpose()
//	}
//	private func minorMatrix(removingRow rowToRemove: Int, removingColumn columnToRemove: Int) -> [[T]] {
//		values.enumerated().compactMap { i, row in
//			guard i != rowToRemove else { return nil }
//			let filteredRow = row.enumerated().compactMap { j, value in
//				j != columnToRemove ? value : nil
//			}
//			return filteredRow
//		}
//	}
//	
//	public func removing(rows rowsToRemove: [Int], columns columnsToRemove: [Int]) throws -> Matrix<T> {
//		for row in rowsToRemove {
//			if row < 0 || row >= rows { throw MatrixError.indexOutOfRange }
//		}
//		for column in columnsToRemove {
//			if column < 0 || column >= columns { throw MatrixError.indexOutOfRange }
//		}
//		let rowSet = Set(rowsToRemove)
//		let columnSet = Set(columnsToRemove)
//		let filteredValues: [[T]] = values.enumerated().compactMap { (i, row) in
//			guard !rowSet.contains(i) else { return nil }
//			let filteredRow = row.enumerated().compactMap { (j, value) in
//				columnSet.contains(j) ? nil : value
//			}
//			return filteredRow
//		}
//		return Matrix(filteredValues)
//	}
//	
//	/// Recursive helper function to calculate the determinant.
//	private static func calculateDeterminant(for matrix: [[T]]) -> T {
//		let n = matrix.count
//	
//		if n == 1 { // 1x1 matrix
//			return matrix[0][0]
//		}
//		if n == 2 { // 2x2 matrix
//			return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]
//		}
//		if n == 3 { // 3x3 matrix
//			let a = matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1])
//			let b = matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0])
//			let c = matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0])
//			return a - b + c
//		}
//		if n < 6 {
//			// Calculate determinant using cofactor expansion - O(n!)
//			var det: T = 0
//			for i in 0..<n {
//				// Create the submatrix by excluding the first row and the ith column
//				var subMatrix = [[T]]()
//				for j in 1..<matrix.count {
//					var row = [T]()
//					for k in 0..<matrix.count {
//						if k != i {
//							row.append(matrix[j][k])
//						}
//					}
//					subMatrix.append(row)
//				}
//				// Recursive determinant calculation with cofactor expansion
//				let sign: T = (i % 2 == 0) ? 1 : -1
//				det += sign * matrix[0][i] * calculateDeterminant(for: subMatrix)
//			}
//			return det
//		} else {
//			// Calculate determinant using LU decomposition - O(n^3)
//			do {
//				let (L, U, _, rowSwaps) = try Matrix.luDecomposition(of: matrix)
//				return calculateDeterminant(fromL: L, withU: U, rowSwaps: rowSwaps)
//			} catch {
//				if let matrixError = error as? MatrixError {
//					if matrixError == MatrixError.singularMatrix { return 0 }
//				}
//				fatalError(error.localizedDescription)
//			}
//		}
//	}
//	/// Calculate the determinant from LU Decomposition
//	private static func calculateDeterminant(fromL L: [[T]], withU U: [[T]], rowSwaps: Int) -> T {
//		//let detL: T = 1
//		var detU: T = 1
//		// Product of diagonal elements
//		for i in 0..<L.count {
//			//detL *= L[i][i] // Diagonal elements are all 1.
//			detU *= U[i][i]
//		}
//		//let det = detL*detU
//		// Adjust sign for number of row swaps
//		return (rowSwaps % 2 == 0) ? detU : -detU
//	}
//	
//	public static func solve<let n: Int>(A: InlineMatrix<n, n, bufferSize, T>, b: InlineMatrix<n, 1, n, T>) throws -> InlineMatrix<n, 1, n, T> {
//		if n < 6 {
//			let inverseA = try A.inverse()
//			return inverseA.multiplying(by: b)
//		}
//		let (L, U, P, _) = try InlineMatrix.luDecomposition(of: A.values)
//		let x = try Matrix.solve(L: L, U: U, P: P, b: Array(b.flatValues))
//		return InlineMatrix<n,1,n,T>(x)
//	}
//	public static func solve(L: [[T]], U: [[T]], P: [Int], b: [T]) throws -> [T] {
//		guard L.count == L.first?.count && U.count == U.first?.count else {
//			throw MatrixError.notSquareMatrix
//		}
//		guard b.count == L.count && L.count == U.count && P.count == L.count else {
//			throw MatrixError.nonMatchingDimensions
//		}
//		let n = L.count
//		
//		// Helper to apply permutation to identity column
//		func permute(_ e: [T]) -> [T] {
//			var result = Array<T>(repeating: 0.0, count: n)
//			for i in 0..<n {
//				result[i] = e[P[i]]
//			}
//			return result
//		}
//		
//		let Pb = permute(b)
//		
//		// Forward substitution: solve L * y = Pb
//		var y = Array<T>(repeating: 0.0, count: n)
//		for j in 0..<n {
//			y[j] = Pb[j]
//			for k in 0..<j {
//				y[j] -= L[j][k] * y[k]
//			}
//			// L diagonal is 1, so divide unnecessary
//		}
//		
//		// Backward substitution: solve U * x = y
//		var x = Array<T>(repeating: 0.0, count: n)
//		for j in stride(from: n - 1, through: 0, by: -1) {
//			x[j] = y[j]
//			for k in (j + 1)..<n {
//				x[j] -= U[j][k] * x[k]
//			}
//			guard U[j][j] != 0 else { throw MatrixError.singularMatrix } // Singular matrix
//			x[j] /= U[j][j]
//		}
//		return x
//	}
//	
//	public func inverse() throws -> InlineMatrix<rows, columns, bufferSize, T> where rows == columns {
//		let n = rows
//		if n < 6 {
//			let det = InlineMatrix.calculateDeterminant(for: values)
//			if det == 0 { throw MatrixError.singularMatrix }
//			
//			if n == 1 { return InlineMatrix([[(1/self[0,0])]]) }
//			if n == 3 {
//				let m = values
//
//				let a = m[0][0], b = m[0][1], c = m[0][2]
//				let d = m[1][0], e = m[1][1], f = m[1][2]
//				let g = m[2][0], h = m[2][1], i = m[2][2]
//
//				let A = e * i - f * h
//				let B = -(d * i - f * g)
//				let C = d * h - e * g
//				let D = -(b * i - c * h)
//				let E = a * i - c * g
//				let F = -(a * h - b * g)
//				let G = b * f - c * e
//				let H = -(a * f - c * d)
//				let I = a * e - b * d
//				
//				//let det = a * A + b * B + c * C
//				let resultValues: [[T]] = [
//					[A, D, G],
//					[B, E, H],
//					[C, F, I]
//				]
//				return InlineMatrix(resultValues)/det
//			}
//			// Use Cofactor Expansion for small matrices (n! < n^3).
//			return (1/det)*(try! self.adjoint())
//		}
//		
//		// TODO: Implement Cholesky decomposition for symmetric positive definite matrices
//		// Use LU Decomposition
//		let (L, U, P, _) = try InlineMatrix.luDecomposition(of: self.values)
//		var invCols: [[T]] = []
//		
//		// Helper to apply permutation to identity column
//		func permute(_ e: [T]) -> [T] {
//			var result = Array<T>(repeating: 0.0, count: n)
//			for i in 0..<n {
//				result[i] = e[P[i]]
//			}
//			return result
//		}
//		
//		// Solve for each column of the inverse
//		for i in 0..<n {
//			// Unit vector (column of identity)
//			var e = Array<T>(repeating: 0.0, count: n)
//			e[i] = 1.0
//			let b = permute(e)
//			
//			// Forward substitution: solve L * y = b
//			var y = Array<T>(repeating: 0.0, count: n)
//			for j in 0..<n {
//				y[j] = b[j]
//				for k in 0..<j {
//					y[j] -= L[j][k] * y[k]
//				}
//				// L diagonal is 1, so divide unnecessary
//			}
//			
//			// Backward substitution: solve U * x = y
//			var x = Array<T>(repeating: 0.0, count: n)
//			for j in stride(from: n - 1, through: 0, by: -1) {
//				x[j] = y[j]
//				for k in (j + 1)..<n {
//					x[j] -= U[j][k] * x[k]
//				}
//				guard U[j][j] != 0 else { throw MatrixError.singularMatrix } // Singular matrix
//				x[j] /= U[j][j]
//			}
//			
//			invCols.append(x)
//		}
//		
//		// Transpose columns to rows
//		let inverseValues = (0..<n).map { i in
//			invCols.map { $0[i] }
//		}
//		return InlineMatrix(inverseValues)
//	}
//	
//	/// Returns (L, U, P, swapCount) for LU decomposition with partial pivoting.
//	/// - Returns: Optional tuple (L, U, P, swapCount), or nil if decomposition fails.
//	public static func luDecomposition(of matrix: [[T]], tolerance: T = 0.000001) throws -> (L: [[T]], U: [[T]], P: [Int], swapCount: Int) {
//		guard matrix.count == matrix.first?.count else { throw MatrixError.notSquareMatrix }
//		let n = matrix.count
//		var A = matrix
//		var L = Array(repeating: Array<T>(repeating: 0.0, count: n), count: n)
//		var U = Array(repeating: Array<T>(repeating: 0.0, count: n), count: n)
//		var P: Array<Int> = Array(0..<n)
//		var swapCount = 0
//		
//		let swapsAreRequired: Bool = {
//			for i in 0..<n {
//				if abs(A[i][i]) < tolerance { return true }
//			}
//			return false
//		}()
//		
//		if swapsAreRequired {
//			for i in 0..<n {
//				// Partial Pivoting
//				var maxRow = i
//				var maxVal: T = 0.0
//				for k in i..<n {
//					let kABS = abs(A[k][i])
//					if kABS > maxVal {
//						maxVal = kABS
//						maxRow = k
//					}
//				}
//				
//				// Zero pivot => singular (or near singular) matrix
//				if maxVal < tolerance {
//					throw MatrixError.singularMatrix
//				}
//				
//				// Swap rows in A and P
//				if maxRow != i {
//					A.swapAt(i, maxRow) // Pivot A
//					P.swapAt(i, maxRow) // Pivot P
//					swapCount += 1
//				}
//			}
//		}
//
//		// Doolittle algorithm + Partial Pivoting
//		for i in 0..<n {
//			
//			// Upper Triangular
//			for j in i..<n {
//				// Summation of L(i, k) * U(k, j)
//				var sum: T = 0.0
//				for k in 0..<i {
//					sum += L[i][k] * U[k][j]
//				}
//				// Evaluating U(i, j)
//				U[i][j] = A[i][j] - sum
//			}
//
//			
//			// Lower Triangular
//			L[i][i] = 1.0 // unit diagonal
//			for j in (i+1)..<n {
//				// Summation of L(j, k) * U(k, i)
//				var sum: T = 0.0
//				for k in 0..<i {
//					sum += L[j][k] * U[k][i]
//				}
//				// Evaluating L(j, i)
//				if U[i][i] == 0.0 { throw MatrixError.singularMatrix } // Singular matrix
//				L[j][i] = (A[j][i] - sum) / U[i][i]
//			}
//		}
//		return (L, U, P, swapCount)
//	}
//	/// Builds a permutation matrix from a row permutation array.
//	/// Each index `i` in `permutation` maps original row `i` to `permutation[i]`.
//	public static func permutationMatrix(from permutation: [Int]) -> InlineMatrix<rows, columns, bufferSize, T> {
//		let n = permutation.count
//		var values = Array(repeating: Array<T>(repeating: 0.0, count: n), count: n)
//		for i in 0..<n {
//			values[permutation[i]][i] = 1.0
//		}
//		return InlineMatrix(values)
//	}
//	
//	// UNTESTED - written by AI - LU decomposition is about 2x faster than QR usually for solving matrices, but QR is best for eigenvalues and eigenvectors.
//	private func qrDecompose(iterations: Int = 10) -> (Q: [[T]], R: [[T]]) {
//		var A = self.values
//		let n = A.count
//		let m = A[0].count
//		var Q: [[T]] = Array(repeating: Array(repeating: 0.0, count: m), count: n)
//		var R: [[T]] = Array(repeating: Array(repeating: 0.0, count: m), count: m)
//
//		for k in 0..<m {
//			var norm: T = 0.0
//			for i in 0..<n {
//				norm += A[i][k] * A[i][k]
//			}
//			norm = sqrt(norm)
//			for i in 0..<n {
//				Q[i][k] = A[i][k] / norm
//			}
//			R[k][k] = norm
//
//			for j in (k+1)..<m {
//				var dot: T = 0.0
//				for i in 0..<n {
//					dot += Q[i][k] * A[i][j]
//				}
//				R[k][j] = dot
//				for i in 0..<n {
//					A[i][j] -= Q[i][k] * dot
//				}
//			}
//		}
//		return (Q, R)
//	}
//	
//	/// Display row count x column count
//	public var sizeDescription: String {
//		"\(rows)x\(columns)"
//	}
//	
//	/// Convert matrix to a string representation
//	public var description: String {
//		var result = "\(rows)x\(columns)\n"
//		for i in 0..<rows {
//			result += "["
//			for j in 0..<columns {
//				result += "\(self[i,j])"
//				if j < columns - 1 {
//					result += ", "
//				}
//			}
//			result += "]\n"
//		}
//		return result
//	}
//	
//	/// Check if indices are valid
//	private func isValidIndex(row: Int, col: Int) -> Bool {
//		return row >= 0 && row < rows && col >= 0 && col < columns
//	}
//	/// Returns true if input dimensions are valid.
//	public func isValidMatrix() -> Bool {
//		rows*columns == bufferSize && rows > 0 && columns > 0
//	}
////	/// Fetches or updates a value in the matrix.
////	public subscript(safely r: Int, c: Int) -> T? {
////		get {
////			guard isValidIndex(row: r, col: c) else {
////				return nil
////				print("Matrix index out of bounds: (\(r), \(c))")
////			}
////			return flatValues[unchecked: r * columns + c]
////		}
////		set(newValue) {
////			guard isValidIndex(row: r, col: c) else {
////				print("Matrix index out of bounds: (\(r), \(c))")
////			}
////			if let newValue {
////				flatValues[unchecked: r * columns + c] = newValue
////			}
////		}
////	}
//	/// Unsafely fetches or updates a value in the matrix.
//	private subscript(r: Int, c: Int) -> T {
//		get {
//			return flatValues[unchecked: r * columns + c]
//		}
//		set(newValue) {
//			flatValues[unchecked: r * columns + c] = newValue
//		}
//	}
//	
//	// MARK: Surge
//	// Surge MIT License:
//
//	//Copyright © 2014-2019 the Surge contributors
//	//
//	//Permission is hereby granted, free of charge, to any person obtaining a copy
//	//of this software and associated documentation files (the "Software"), to deal
//	//in the Software without restriction, including without limitation the rights
//	//to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
//	//copies of the Software, and to permit persons to whom the Software is
//	//furnished to do so, subject to the following conditions:
//	//
//	//The above copyright notice and this permission notice shall be included in
//	//all copies or substantial portions of the Software.
//	//
//	//THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//	//IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//	//FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
//	//AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//	//LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
//	//OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
//	//THE SOFTWARE.
//	
//	public func eigenvalues(maxIterations: Int = 100, tolerance: Double = 1e-10) throws -> [(x: T, i: T)] {
//		guard isSquare else { throw MatrixError.notSquareMatrix }
//		if rows == 1 {
//			return self[0,0].isNaN ? [] : [((self[0,0]), T(0))]
//		}
//		if rows == 2 {
//			let a = self[0,0]
//			let b = self[0,1]
//			let c = self[1,0]
//			let d = self[1,1]
//
//			let trace = a + d
//			let determinant = a * d - b * c
//			let discriminant = trace * trace - 4 * determinant
//
//			if discriminant >= 0 {
//				let sqrtDiscriminant = discriminant.squareRoot()
//				let lambda1 = (trace + sqrtDiscriminant) / 2
//				let lambda2 = (trace - sqrtDiscriminant) / 2
//				return [(lambda1,0.0), (lambda2,0.0)]
//			}
//		}
//		return try eigenDecompose(computeEigenVectors: false).eigenValues
//	}
//	
//	/// Decomposes a square matrix into its eigenvalues and left and right eigenvectors.
//	/// The decomposition may result in complex numbers, represented by (T, T), which
//	///   are the (real, imaginary) parts of the complex number.
//	/// - Returns: a struct with the eigen values and left and right eigen vectors using (T, T)
//	///   to represent a complex number.
//	func eigenDecompose(computeEigenVectors: Bool) throws -> MatrixEigenDecompositionResult<T> {
//		let input = Matrix<Double>(self)
//		let decomposition = try eigenDecomposeOLD(input, computeEigenVectors: computeEigenVectors)
//		
////		guard rows == columns else {
////			throw MatrixError.notSquareMatrix
////		}
////		let (Q,R) = self.qrDecompose()
////		let lambdas: [(x: T, i: T)] = []
////		var rightEigenVectors: [[(T, T)]] = [[]]
////		if computeEigenVectors  {
////			for lambda in lambdas {
////				var A = try self-Matrix<T>.identity(size: rows)*lambda.x
////				let realVector = try Matrix.solve(A: A, b: Matrix(Array<T>(repeating: 0, count: rows), isRow: false))
////				A = try self-Matrix<T>.identity(size: rows)*lambda.i
////				let imaginaryVector = try Matrix.solve(A: A, b: Matrix(Array<T>(repeating: 0, count: rows), isRow: false))
////				var eigenVector: [(T, T)] = []
////				for i in 0..<rows {
////					eigenVector.append((realVector.flatValues[i],imaginaryVector.flatValues[i]))
////				}
////				rightEigenVectors.append(eigenVector)
////			}
////		}
//		
//		return MatrixEigenDecompositionResult<T>(
//			eigenValues: decomposition.eigenValues.map { (T($0.0), T($0.1)) },
//			leftEigenVectors: decomposition.leftEigenVectors.map { $0.map { (T($0.0), T($0.1)) } },
//			rightEigenVectors: decomposition.rightEigenVectors.map { $0.map { (T($0.0), T($0.1)) } }
//		)
//	}
//
//	/// Decomposes a square matrix into its eigenvalues and left and right eigenvectors.
//	/// The decomposition may result in complex numbers, represented by (Double, Double), which
//	///   are the (real, imaginary) parts of the complex number.
//	/// - Parameters:
//	///   - lhs: a square matrix
//	/// - Returns: a struct with the eigen values and left and right eigen vectors using (Double, Double)
//	///   to represent a complex number.
//	func eigenDecomposeOLD(_ lhs: Matrix<Double>, computeEigenVectors: Bool) throws -> MatrixEigenDecompositionResult<Double> {
//		
//		guard lhs.rows == lhs.columns else {
//			throw MatrixError.notSquareMatrix
//		}
//
//		// dgeev_ needs column-major matrices, so transpose 'lhs'.
//		var matrixGrid: [__CLPK_doublereal] = lhs.transpose().flatValues
//		var matrixRowCount = __CLPK_integer(lhs.rows)
//		let matrixColCount = matrixRowCount
//		var eigenValueCount = matrixRowCount
//		var leftEigenVectorCount = matrixRowCount
//		var rightEigenVectorCount = matrixRowCount
//
//		var workspaceQuery: Double = 0.0
//		var workspaceSize = __CLPK_integer(-1)
//		var error: __CLPK_integer = 0
//
//		var eigenValueRealParts = [Double](repeating: 0, count: Int(eigenValueCount))
//		var eigenValueImaginaryParts = [Double](repeating: 0, count: Int(eigenValueCount))
//		var leftEigenVectorWork = [Double](repeating: 0, count: Int(leftEigenVectorCount * matrixColCount))
//		var rightEigenVectorWork = [Double](repeating: 0, count: Int(rightEigenVectorCount * matrixColCount))
//
//		var decompositionJobVL: [CChar]
//		var decompositionJobVR: [CChar]
//		
//		if computeEigenVectors  {
//			decompositionJobVL = [0x56, 0x00] // "V" (compute)
//			decompositionJobVR = [0x56, 0x00] // "V" (compute)
//		} else {
//			decompositionJobVL = Array("N".utf8CString) // "N" (do not compute)
//			decompositionJobVR = Array("N".utf8CString) // "N" (do not compute)
//		}
//
//		// Call dgeev to find out how much workspace to allocate
//		dgeev_(&decompositionJobVL, &decompositionJobVR, &matrixRowCount, &matrixGrid, &eigenValueCount, &eigenValueRealParts, &eigenValueImaginaryParts, &leftEigenVectorWork, &leftEigenVectorCount, &rightEigenVectorWork, &rightEigenVectorCount, &workspaceQuery, &workspaceSize, &error)
//		if error != 0 {
//			throw MatrixError.matrixNotDecomposable
//		}
//
//		// Allocate the workspace and call dgeev again to do the actual decomposition
//		var workspace = [Double](repeating: 0.0, count: Int(workspaceQuery))
//		workspaceSize = __CLPK_integer(workspaceQuery)
//		dgeev_(&decompositionJobVL, &decompositionJobVR, &matrixRowCount, &matrixGrid, &eigenValueCount, &eigenValueRealParts, &eigenValueImaginaryParts, &leftEigenVectorWork, &leftEigenVectorCount, &rightEigenVectorWork, &rightEigenVectorCount, &workspace, &workspaceSize, &error)
//		if error != 0 {
//			throw MatrixError.matrixNotDecomposable
//		}
//
//		return MatrixEigenDecompositionResult<Double>(rowCount: lhs.rows, eigenValueRealParts: eigenValueRealParts, eigenValueImaginaryParts: eigenValueImaginaryParts, leftEigenVectorWork: leftEigenVectorWork, rightEigenVectorWork: rightEigenVectorWork)
//	}
//	
//	/// Holds the result of eigendecomposition. The (Scalar, Scalar) used
//	/// in the property types represents a complex number with (real, imaginary) parts.
//	struct MatrixEigenDecompositionResult<Scalar: BinaryFloatingPoint> {
//		public let eigenValues: [(Scalar, Scalar)]
//		public let leftEigenVectors: [[(Scalar, Scalar)]]
//		public let rightEigenVectors: [[(Scalar, Scalar)]]
//		
//		public init(eigenValues: [(Scalar, Scalar)], leftEigenVectors: [[(Scalar, Scalar)]], rightEigenVectors: [[(Scalar, Scalar)]]) {
//			self.eigenValues = eigenValues
//			self.leftEigenVectors = leftEigenVectors
//			self.rightEigenVectors = rightEigenVectors
//		}
//		
//		public init(rowCount: Int, eigenValueRealParts: [Scalar], eigenValueImaginaryParts: [Scalar], leftEigenVectorWork: [Scalar], rightEigenVectorWork: [Scalar]) {
//			// The eigenvalues are an array of (real, imaginary) results from dgeev
//			self.eigenValues = Array(zip(eigenValueRealParts, eigenValueImaginaryParts))
//			
//			// Build the left and right eigenvectors
//			let emptyVector = [(Scalar, Scalar)](repeating: (0.0, 0.0), count: rowCount)
//			var leftEigenVectors = [[(Scalar, Scalar)]](repeating: emptyVector, count: rowCount)
//			MatrixEigenDecompositionResult.buildEigenVector(eigenValueImaginaryParts: eigenValueImaginaryParts, eigenVectorWork: leftEigenVectorWork, result: &leftEigenVectors)
//			
//			var rightEigenVectors = [[(Scalar, Scalar)]](repeating: emptyVector, count: rowCount)
//			MatrixEigenDecompositionResult.buildEigenVector(eigenValueImaginaryParts: eigenValueImaginaryParts, eigenVectorWork: rightEigenVectorWork, result: &rightEigenVectors)
//			
//			self.leftEigenVectors = leftEigenVectors
//			self.rightEigenVectors = rightEigenVectors
//		}
//		
//		// Convert the result of dgeev into an array of complex numbers
//		// See Intel's documentation on column-major results for sample code that this is based on:
//		// https://software.intel.com/sites/products/documentation/doclib/mkl_sa/11/mkl_lapack_examples/dgeev.htm
//		private static func buildEigenVector<S>(eigenValueImaginaryParts: [S], eigenVectorWork: [S], result: inout [[(S, S)]]) where S: FloatingPoint & ExpressibleByFloatLiteral {
//			// row and col count are the same because result must be square.
//			let rowColCount = result.count
//
//			for row in 0..<rowColCount {
//				var col = 0
//				while col < rowColCount {
//					if eigenValueImaginaryParts[col] == 0.0 {
//						// v is column-major
//						result[row][col] = (eigenVectorWork[row + rowColCount * col], 0.0)
//						col += 1
//					} else {
//						// v is column-major
//						result[row][col] = (eigenVectorWork[row + col * rowColCount], eigenVectorWork[row + rowColCount * (col + 1)])
//						result[row][col + 1] = (eigenVectorWork[row + col * rowColCount], -eigenVectorWork[row + rowColCount * (col + 1)])
//						col += 2
//					}
//				}
//			}
//		}
//
//	}
//}

//@available(iOS 26.0.0, macOS 26.0.0, *)
//nonisolated private struct InlineSquareMatrix<let n: Int, let bufferSize: Int, T: BinaryFloatingPoint>: Equatable, CustomStringConvertible, Copyable {
//	
//	public var values: [[T]] {
//		get {
//			var result = [[T]]()
//			for i in 0..<bufferSize {
//				let r = i / n
//				let c = i % n
//				result[r][c] = flatValues[i]
//			}
//			return result
//		}
//		set {
//			for i in 0..<bufferSize {
//				let r = i / n
//				let c = i % n
//				flatValues[i] = newValue[r][c]
//			}
//		}
//	}
//	public private(set)var flatValues: InlineArray<bufferSize, T>
//	public var rows: Int { n }
//	public var columns: Int { n }
//	
//	static func == (lhs: InlineSquareMatrix<n, bufferSize, T>, rhs: InlineSquareMatrix<n, bufferSize, T>) -> Bool {
//		lhs.rows == rhs.rows && lhs.flatValues == rhs.flatValues
//	}
//	func isIdentical(to rhs: InlineSquareMatrix<n, bufferSize, T>) -> Bool {
//		self.flatValues.span.isIdentical(to: rhs.flatValues.span)
//	}
//	
//	public init<let r: Int, let c: Int, V>(_ matrix: InlineMatrix<r, c, bufferSize, V>) where V: BinaryFloatingPoint {
//		self.flatValues = InlineArray<bufferSize, T>() { T(matrix.flatValues[$0]) }
//	}
//	public init<V>(_ matrix: InlineSquareMatrix<n, bufferSize, V>) where V: BinaryFloatingPoint {
//		self.flatValues = InlineArray<bufferSize, T>() { T(matrix.flatValues[$0]) }
//	}
//	public init(_ matrix: InlineSquareMatrix<n, bufferSize, T>) {
//		self.flatValues = matrix.flatValues
//	}
//	/// Creates a matrix where the outer array is for rows and the inner array is for columns.
//	public init(_ values: [[T]]) {
//		let tempFlatValues = values.flatMap(\.self)
//		self.flatValues = InlineArray() {
//			tempFlatValues[$0]
//		}
//	}
//	public init(_ values: [T]) {
//		self.flatValues = InlineArray() { values[$0] }
//	}
//	public init(_ values: InlineArray<bufferSize, T>) {
//		self.flatValues = values
//	}
//	
//	/// Initialize a matrix with a default value
//	public init(repeating value: T = 0) {
//		self.flatValues = InlineArray<bufferSize, T>(repeating: value)
//	}
//	
//	public init(_ closure: (_ row: Int, _ column: Int) throws -> T) rethrows {
//		self.flatValues = try InlineArray<bufferSize, T>() { i in
//			let r = i / n
//			let c = i % n
//			return try closure(r, c)
//		}
//	}
//	public init(_ closure: (_ i: Int) throws -> T) rethrows {
//		self.flatValues = try InlineArray<bufferSize, T>() { i in
//			return try closure(i)
//		}
//	}
//
//	public static func identity() -> InlineSquareMatrix<n, bufferSize, T> {
//		return InlineSquareMatrix<n, bufferSize, T>() { i in
//			if i / n == i % n { 1 } else { 0 }
//		}
//	}
//	public static func diagonal(repeating value: T) -> InlineSquareMatrix<n, bufferSize, T> {
//		return InlineSquareMatrix<n, bufferSize, T>() { i in
//			if i / n == i % n { value } else { 0 }
//		}
//	}
//	public static func diagonal<C>(scalars: C) -> InlineSquareMatrix<n, bufferSize, T> where C: Collection, C.Element == T {
//		var iScalar: C.Index = scalars.index(scalars.startIndex, offsetBy: -1)
//		return InlineSquareMatrix<n, bufferSize, T>() { i in
//			if i / n == i % n {
//				iScalar = scalars.index(after: iScalar)
//				return scalars[iScalar]
//			} else { return 0 }
//		}
//	}
//
//	public init(_ simdMatrix: simd_double2x2) where T == Double, n == 2 {
//		let m = simdMatrix.transpose
//		let temp = [
//			m.columns.0.x, m.columns.0.y,
//			m.columns.1.x, m.columns.1.y,
//		]
//		self.flatValues = InlineArray() { temp[$0] }
//	}
//	public init(_ simdMatrix: simd_double3x3) where T == Double, n == 3 {
//		let m = simdMatrix.transpose
//		let temp = [
//			m.columns.0.x, m.columns.0.y, m.columns.0.z,
//			m.columns.1.x, m.columns.1.y, m.columns.1.z,
//			m.columns.2.x, m.columns.2.y, m.columns.2.z,
//		]
//		self.flatValues = InlineArray() { temp[$0] }
//	}
//	public init(_ simdMatrix: simd_double4x4) where T == Double, n == 4{
//		let m = simdMatrix.transpose
//		let temp = [
//			m.columns.0.x, m.columns.0.y, m.columns.0.z, m.columns.0.w,
//			m.columns.1.x, m.columns.1.y, m.columns.1.z, m.columns.1.w,
//			m.columns.2.x, m.columns.2.y, m.columns.2.z, m.columns.2.w,
//			m.columns.3.x, m.columns.3.y, m.columns.3.z, m.columns.3.w,
//		]
//		self.flatValues = InlineArray() { temp[$0] }
//	}
//	
//	public subscript(row row: Int) -> InlineArray<n, T> {
//		get {
//			assert(row < n)
//			let startIndex = row * n
//			return InlineArray<n, T>() { c in
//				flatValues[startIndex+c]
//			}
//		}
//		set {
//			assert(row < n)
//			assert(newValue.count == n)
//			let startIndex = row * n
//			for c in 0..<newValue.count {
//				flatValues[c+startIndex] = newValue[c]
//			}
//		}
//	}
//	public subscript(column column: Int) -> InlineArray<n, T> {
//		get {
//			return InlineArray<n, T>() { r in
//				let index = r * n + column
//				return self.flatValues[index]
//			}
//		}
//		set {
//			assert(column < n)
//			assert(newValue.count == n)
//			for r in 0..<n {
//				let index = r * n + column
//				flatValues[index] = newValue[r]
//			}
//		}
//	}
//	
//	/// Returns `true` if both matrices are approximately equivalent to 10^-4
//	public static func ~= (lhs: InlineSquareMatrix<n, bufferSize, T>, rhs: InlineSquareMatrix<n, bufferSize,T>) -> Bool {
//		let result = lhs - rhs
//		for i in 0..<bufferSize {
//			if abs(result.flatValues[i]) > 0.000001 {
//				return false
//			}
//		}
//		return true
//	}
//	
//	/// Elementwise addition
//	public static func + (lhs: InlineSquareMatrix<n, bufferSize,T>, rhs: InlineSquareMatrix<n, bufferSize,T>) -> InlineSquareMatrix<n, bufferSize,T> {
//		var result = lhs
//		for i in 0..<bufferSize {
//			result.flatValues[i] = lhs.flatValues[i] - rhs.flatValues[i]
//		}
//		return result
//	}
//	/// Elementwise subtraction
//	public static func - (lhs: InlineSquareMatrix<n, bufferSize,T>, rhs: InlineSquareMatrix<n, bufferSize,T>) -> InlineSquareMatrix<n, bufferSize,T> {
//		var result = lhs
//		for i in 0..<bufferSize {
//			result.flatValues[i] = lhs.flatValues[i] - rhs.flatValues[i]
//		}
//		return result
//	}
//	/// Scalar multiplication
//	public static func * (lhs: InlineSquareMatrix<n, bufferSize,T>, rhs: T) -> InlineSquareMatrix<n, bufferSize,T> {
//		var result = lhs
//		for i in 0..<bufferSize {
//			result.flatValues[i] = lhs.flatValues[i] * rhs
//		}
//		return result
//	}
//	/// Scalar division
//	public static func / (lhs: InlineSquareMatrix<n, bufferSize,T>, rhs: T) -> InlineSquareMatrix<n, bufferSize,T> {
//		var result = lhs
//		for i in 0..<bufferSize {
//			result.flatValues[i] = lhs.flatValues[i] / rhs
//		}
//		return result
//	}
//	/// Scalar multiplication
//	public static func * (lhs: T, rhs: InlineSquareMatrix<n, bufferSize,T>) -> InlineSquareMatrix<n, bufferSize,T> {
//		var result = rhs
//		for i in 0..<bufferSize {
//			result.flatValues[i] = rhs.flatValues[i] * lhs
//		}
//		return result
//	}
//	/// Elementwise multiplication of two matricies
//	public mutating func multiplyElementwise(by other: Self) {
//		for i in 0..<bufferSize {
//			flatValues[i] = flatValues[i]*other.flatValues[i]
//		}
//	}
//	/// Elementwise multiplication of two matricies
//	public func multiplyingElementwise(by other: Self) -> Self {
//		var newMatrix = self
//		newMatrix.multiplyElementwise(by: other)
//		return newMatrix
//	}
////	/// Perform Matrix Multiplication
////	public func multiplying<let M: Int, let S: Int, let S2: Int>(by rhs: InlineMatrix<n,M,S,T>) -> InlineMatrix<n, M, S2, T> {
////		var result = InlineMatrix<n, M, S2, T>()
////		for i in 0..<n { // lhs rows
////			for j in 0..<M { // rhs columns
////				for k in 0..<n { // lhs colmuns
////					result[i,j] += self[i,k] * rhs[k,j]
////				}
////			}
////		}
////		return result
////	}
//	/// Perform Matrix Multiplication
//	public func multiplying(by rhs: InlineArray<n,T>) -> InlineArray<n, T> {
//		var result = InlineArray<n, T>(repeating: T.zero)
//		for i in 0..<n { // lhs rows
//			for k in 0..<n { // lhs colmuns
//				result[i] += self[i,k] * rhs[k]
//			}
//		}
//		return result
//	}
//	
//	/// Returns the matrix as a single inline array, only if matrix is a single column or a single row.
//	public func makeInlineVector() throws -> InlineArray<bufferSize,T> {
//		guard rows == 1 || columns == 1 else { throw MatrixError.nonMatchingDimensions }
//		if rows == 1 {
//			guard columns == bufferSize else { throw MatrixError.nonMatchingDimensions }
//		}
//		if columns == 1 {
//			guard rows == bufferSize else { throw MatrixError.nonMatchingDimensions }
//		}
//		return flatValues
//	}
//	/// Returns the matrix as a single array, only if matrix is a single column or a single row.
//	public func makeVectorArray() throws -> Array<T> {
//		return try Array(makeInlineVector())
//	}
//	
//	/// Returns the transpose of the matrix (swapping rows for columns)
//	public func transpose() -> InlineSquareMatrix<n, bufferSize, T> {
//		return InlineSquareMatrix<n, bufferSize, T> { row, column in
//			self[column, row]
//		}
//	}
//	
//	/// Returns `true` if matrix is symmetric
//	public var isSymmetric: Bool {
//		for r in 0..<n {
//			for c in 0..<n {
//				let i1 = r * n + c
//				let i2 = c * n + r
//				if flatValues[i1] != flatValues[i2] { return false }
//			}
//		}
//		return true
//	}
//	
//	/// Determinant of the matrix (only defined for square matrices).
//	public func determinant() throws -> T {
//		if n == 1 { // 1x1 matrix
//			return self[0,0]
//		}
//		if n == 2 { // 2x2 matrix
//			return self[0,0] * self[1,1] - self[0,1] * self[1,0]
//		}
//		if n == 3 { // 3x3 matrix
//			let a = self[0,0] * (self[1,1] * self[2,2] - self[1,2] * self[2,1])
//			let b = self[0,1] * (self[1,0] * self[2,2] - self[1,2] * self[2,0])
//			let c = self[0,2] * (self[1,0] * self[2,1] - self[1,1] * self[2,0])
//			return a - b + c
//		}
//		if n < 6 {
//			// Calculate determinant using cofactor expansion - O(n!)
//			var det: T = 0
//			for i in 0..<n {
//				// Create the submatrix by excluding the first row and the ith column
//				var subMatrix = [[T]]()
//				for j in 1..<n {
//					var row = [T]()
//					for k in 0..<n {
//						if k != i {
//							row.append(self[j,k])
//						}
//					}
//					subMatrix.append(row)
//				}
//				// Recursive determinant calculation with cofactor expansion
//				let sign: T = (i % 2 == 0) ? 1 : -1
//				det += sign * self[0,i] * InlineSquareMatrix.calculateDeterminant(for: subMatrix)
//			}
//			return det
//		} else {
//			// Calculate determinant using LU decomposition - O(n^3)
//			do {
//				let (L, U, _, rowSwaps) = try Matrix.luDecomposition(of: values)
//				return InlineSquareMatrix.calculateDeterminant(fromL: L, withU: U, rowSwaps: rowSwaps)
//			} catch {
//				if let matrixError = error as? MatrixError {
//					if matrixError == MatrixError.singularMatrix { return 0 }
//				}
//				fatalError(error.localizedDescription)
//			}
//		}
//	}
//	
//	public func cofactorMatrix() throws -> InlineSquareMatrix<n, bufferSize, T> {
//		guard n > 1 else { throw MatrixError.indexOutOfRange }
//		if n == 2 {
//			return InlineSquareMatrix([
//				self[1,1], -self[1,0],
//				-self[0,1], self[0,0]
//			])
//		}
//		var cofactorMatrix = InlineSquareMatrix<n, bufferSize, T>() { r, c in
//			let minor = minorMatrix(removingRow: r, removingColumn: c)
//			let sign: T = ((r + c) % 2 == 0) ? 1 : -1
//			cofactorMatrix[r,c] = sign * InlineSquareMatrix.calculateDeterminant(for: minor)
//		}
//		return cofactorMatrix
//	}
//	public func adjoint() throws -> InlineSquareMatrix<n, bufferSize, T> {
//		// Transpose the cofactor matrix to get the adjoint
//		return try self.cofactorMatrix().transpose()
//	}
//	private func minorMatrix(removingRow rowToRemove: Int, removingColumn columnToRemove: Int) -> [[T]] {
//		values.enumerated().compactMap { i, row in
//			guard i != rowToRemove else { return nil }
//			let filteredRow = row.enumerated().compactMap { j, value in
//				j != columnToRemove ? value : nil
//			}
//			return filteredRow
//		}
//	}
//	
//	public func removing(rows rowsToRemove: [Int], columns columnsToRemove: [Int]) throws -> Matrix<T> {
//		for row in rowsToRemove {
//			if row < 0 || row >= rows { throw MatrixError.indexOutOfRange }
//		}
//		for column in columnsToRemove {
//			if column < 0 || column >= columns { throw MatrixError.indexOutOfRange }
//		}
//		let rowSet = Set(rowsToRemove)
//		let columnSet = Set(columnsToRemove)
//		let filteredValues: [[T]] = values.enumerated().compactMap { (i, row) in
//			guard !rowSet.contains(i) else { return nil }
//			let filteredRow = row.enumerated().compactMap { (j, value) in
//				columnSet.contains(j) ? nil : value
//			}
//			return filteredRow
//		}
//		return Matrix(filteredValues)
//	}
//	
//	/// Recursive helper function to calculate the determinant.
//	private static func calculateDeterminant(for matrix: [[T]]) -> T {
//		let n = matrix.count
//	
//		if n == 1 { // 1x1 matrix
//			return matrix[0][0]
//		}
//		if n == 2 { // 2x2 matrix
//			return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]
//		}
//		if n == 3 { // 3x3 matrix
//			let a = matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1])
//			let b = matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0])
//			let c = matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0])
//			return a - b + c
//		}
//		if n < 6 {
//			// Calculate determinant using cofactor expansion - O(n!)
//			var det: T = 0
//			for i in 0..<n {
//				// Create the submatrix by excluding the first row and the ith column
//				var subMatrix = [[T]]()
//				for j in 1..<matrix.count {
//					var row = [T]()
//					for k in 0..<matrix.count {
//						if k != i {
//							row.append(matrix[j][k])
//						}
//					}
//					subMatrix.append(row)
//				}
//				// Recursive determinant calculation with cofactor expansion
//				let sign: T = (i % 2 == 0) ? 1 : -1
//				det += sign * matrix[0][i] * calculateDeterminant(for: subMatrix)
//			}
//			return det
//		} else {
//			// Calculate determinant using LU decomposition - O(n^3)
//			do {
//				let (L, U, _, rowSwaps) = try Matrix.luDecomposition(of: matrix)
//				return calculateDeterminant(fromL: L, withU: U, rowSwaps: rowSwaps)
//			} catch {
//				if let matrixError = error as? MatrixError {
//					if matrixError == MatrixError.singularMatrix { return 0 }
//				}
//				fatalError(error.localizedDescription)
//			}
//		}
//	}
//	/// Calculate the determinant from LU Decomposition
//	private static func calculateDeterminant(fromL L: [[T]], withU U: [[T]], rowSwaps: Int) -> T {
//		//let detL: T = 1
//		var detU: T = 1
//		// Product of diagonal elements
//		for i in 0..<L.count {
//			//detL *= L[i][i] // Diagonal elements are all 1.
//			detU *= U[i][i]
//		}
//		//let det = detL*detU
//		// Adjust sign for number of row swaps
//		return (rowSwaps % 2 == 0) ? detU : -detU
//	}
//	
//	public static func solve(A: InlineSquareMatrix<n, bufferSize, T>, b: InlineArray<n, T>) throws -> InlineArray<n, T> {
//		if n < 6 {
//			let inverseA = try A.inverse()
//			return inverseA.multiplying(by: b)
//		}
//		let (L, U, P, _) = try InlineSquareMatrix.luDecomposition(of: A.values)
//		let x = try InlineSquareMatrix.solve(L: L, U: U, P: P, b: Array(b))
//		return InlineArray<n,T>() { x[$0] }
//	}
//	public static func solve(L: [[T]], U: [[T]], P: [Int], b: [T]) throws -> [T] {
//		guard L.count == L.first?.count && U.count == U.first?.count else {
//			throw MatrixError.notSquareMatrix
//		}
//		guard b.count == L.count && L.count == U.count && P.count == L.count else {
//			throw MatrixError.nonMatchingDimensions
//		}
//		let n = L.count
//		
//		// Helper to apply permutation to identity column
//		func permute(_ e: [T]) -> [T] {
//			var result = Array<T>(repeating: 0.0, count: n)
//			for i in 0..<n {
//				result[i] = e[P[i]]
//			}
//			return result
//		}
//		
//		let Pb = permute(b)
//		
//		// Forward substitution: solve L * y = Pb
//		var y = Array<T>(repeating: 0.0, count: n)
//		for j in 0..<n {
//			y[j] = Pb[j]
//			for k in 0..<j {
//				y[j] -= L[j][k] * y[k]
//			}
//			// L diagonal is 1, so divide unnecessary
//		}
//		
//		// Backward substitution: solve U * x = y
//		var x = Array<T>(repeating: 0.0, count: n)
//		for j in stride(from: n - 1, through: 0, by: -1) {
//			x[j] = y[j]
//			for k in (j + 1)..<n {
//				x[j] -= U[j][k] * x[k]
//			}
//			guard U[j][j] != 0 else { throw MatrixError.singularMatrix } // Singular matrix
//			x[j] /= U[j][j]
//		}
//		return x
//	}
//	
//	public func inverse() throws -> InlineSquareMatrix<n, bufferSize, T> {
//		if n < 6 {
//			let det = InlineSquareMatrix.calculateDeterminant(for: values)
//			if det == 0 { throw MatrixError.singularMatrix }
//			
//			if n == 1 { return InlineSquareMatrix([[(1/self[0,0])]]) }
//			if n == 3 {
//				let m = values
//
//				let a = m[0][0], b = m[0][1], c = m[0][2]
//				let d = m[1][0], e = m[1][1], f = m[1][2]
//				let g = m[2][0], h = m[2][1], i = m[2][2]
//
//				let A = e * i - f * h
//				let B = -(d * i - f * g)
//				let C = d * h - e * g
//				let D = -(b * i - c * h)
//				let E = a * i - c * g
//				let F = -(a * h - b * g)
//				let G = b * f - c * e
//				let H = -(a * f - c * d)
//				let I = a * e - b * d
//				
//				//let det = a * A + b * B + c * C
//				let resultValues: [[T]] = [
//					[A, D, G],
//					[B, E, H],
//					[C, F, I]
//				]
//				return InlineSquareMatrix(resultValues)/det
//			}
//			// Use Cofactor Expansion for small matrices (n! < n^3).
//			return (1/det)*(try! self.adjoint())
//		}
//		
//		// TODO: Implement Cholesky decomposition for symmetric positive definite matrices
//		// Use LU Decomposition
//		let (L, U, P, _) = try InlineSquareMatrix.luDecomposition(of: self.values)
//		var invCols: [[T]] = []
//		
//		// Helper to apply permutation to identity column
//		func permute(_ e: [T]) -> [T] {
//			var result = Array<T>(repeating: 0.0, count: n)
//			for i in 0..<n {
//				result[i] = e[P[i]]
//			}
//			return result
//		}
//		
//		// Solve for each column of the inverse
//		for i in 0..<n {
//			// Unit vector (column of identity)
//			var e = Array<T>(repeating: 0.0, count: n)
//			e[i] = 1.0
//			let b = permute(e)
//			
//			// Forward substitution: solve L * y = b
//			var y = Array<T>(repeating: 0.0, count: n)
//			for j in 0..<n {
//				y[j] = b[j]
//				for k in 0..<j {
//					y[j] -= L[j][k] * y[k]
//				}
//				// L diagonal is 1, so divide unnecessary
//			}
//			
//			// Backward substitution: solve U * x = y
//			var x = Array<T>(repeating: 0.0, count: n)
//			for j in stride(from: n - 1, through: 0, by: -1) {
//				x[j] = y[j]
//				for k in (j + 1)..<n {
//					x[j] -= U[j][k] * x[k]
//				}
//				guard U[j][j] != 0 else { throw MatrixError.singularMatrix } // Singular matrix
//				x[j] /= U[j][j]
//			}
//			
//			invCols.append(x)
//		}
//		
//		// Transpose columns to rows
//		let inverseValues = (0..<n).map { i in
//			invCols.map { $0[i] }
//		}
//		return InlineSquareMatrix(inverseValues)
//	}
//	
//	/// Returns (L, U, P, swapCount) for LU decomposition with partial pivoting.
//	/// - Returns: Optional tuple (L, U, P, swapCount), or nil if decomposition fails.
//	public static func luDecomposition(of matrix: [[T]], tolerance: T = 0.000001) throws -> (L: [[T]], U: [[T]], P: [Int], swapCount: Int) {
//		guard matrix.count == matrix.first?.count else { throw MatrixError.notSquareMatrix }
//		let n = matrix.count
//		var A = matrix
//		var L = Array(repeating: Array<T>(repeating: 0.0, count: n), count: n)
//		var U = Array(repeating: Array<T>(repeating: 0.0, count: n), count: n)
//		var P: Array<Int> = Array(0..<n)
//		var swapCount = 0
//		
//		let swapsAreRequired: Bool = {
//			for i in 0..<n {
//				if abs(A[i][i]) < tolerance { return true }
//			}
//			return false
//		}()
//		
//		if swapsAreRequired {
//			for i in 0..<n {
//				// Partial Pivoting
//				var maxRow = i
//				var maxVal: T = 0.0
//				for k in i..<n {
//					let kABS = abs(A[k][i])
//					if kABS > maxVal {
//						maxVal = kABS
//						maxRow = k
//					}
//				}
//				
//				// Zero pivot => singular (or near singular) matrix
//				if maxVal < tolerance {
//					throw MatrixError.singularMatrix
//				}
//				
//				// Swap rows in A and P
//				if maxRow != i {
//					A.swapAt(i, maxRow) // Pivot A
//					P.swapAt(i, maxRow) // Pivot P
//					swapCount += 1
//				}
//			}
//		}
//
//		// Doolittle algorithm + Partial Pivoting
//		for i in 0..<n {
//			
//			// Upper Triangular
//			for j in i..<n {
//				// Summation of L(i, k) * U(k, j)
//				var sum: T = 0.0
//				for k in 0..<i {
//					sum += L[i][k] * U[k][j]
//				}
//				// Evaluating U(i, j)
//				U[i][j] = A[i][j] - sum
//			}
//
//			
//			// Lower Triangular
//			L[i][i] = 1.0 // unit diagonal
//			for j in (i+1)..<n {
//				// Summation of L(j, k) * U(k, i)
//				var sum: T = 0.0
//				for k in 0..<i {
//					sum += L[j][k] * U[k][i]
//				}
//				// Evaluating L(j, i)
//				if U[i][i] == 0.0 { throw MatrixError.singularMatrix } // Singular matrix
//				L[j][i] = (A[j][i] - sum) / U[i][i]
//			}
//		}
//		return (L, U, P, swapCount)
//	}
//	/// Builds a permutation matrix from a row permutation array.
//	/// Each index `i` in `permutation` maps original row `i` to `permutation[i]`.
//	public static func permutationMatrix(from permutation: [Int]) -> InlineSquareMatrix<n, bufferSize, T> {
//		guard n == permutation.count else { fatalError() }
//		var values = Array(repeating: Array<T>(repeating: 0.0, count: n), count: n)
//		for i in 0..<n {
//			values[permutation[i]][i] = 1.0
//		}
//		return InlineSquareMatrix(values)
//	}
//	
//	// UNTESTED - written by AI - LU decomposition is about 2x faster than QR usually for solving matrices, but QR is best for eigenvalues and eigenvectors.
//	private func qrDecompose(iterations: Int = 10) -> (Q: [[T]], R: [[T]]) {
//		var A = self.values
//		let n = A.count
//		let m = A[0].count
//		var Q: [[T]] = Array(repeating: Array(repeating: 0.0, count: m), count: n)
//		var R: [[T]] = Array(repeating: Array(repeating: 0.0, count: m), count: m)
//
//		for k in 0..<m {
//			var norm: T = 0.0
//			for i in 0..<n {
//				norm += A[i][k] * A[i][k]
//			}
//			norm = sqrt(norm)
//			for i in 0..<n {
//				Q[i][k] = A[i][k] / norm
//			}
//			R[k][k] = norm
//
//			for j in (k+1)..<m {
//				var dot: T = 0.0
//				for i in 0..<n {
//					dot += Q[i][k] * A[i][j]
//				}
//				R[k][j] = dot
//				for i in 0..<n {
//					A[i][j] -= Q[i][k] * dot
//				}
//			}
//		}
//		return (Q, R)
//	}
//	
//	/// Display row count x column count
//	public var sizeDescription: String {
//		"\(rows)x\(columns)"
//	}
//	
//	/// Convert matrix to a string representation
//	public var description: String {
//		var result = "\(rows)x\(columns)\n"
//		for i in 0..<rows {
//			result += "["
//			for j in 0..<columns {
//				result += "\(self[i,j])"
//				if j < columns - 1 {
//					result += ", "
//				}
//			}
//			result += "]\n"
//		}
//		return result
//	}
//	
//	/// Check if indices are valid
//	private func isValidIndex(row: Int, col: Int) -> Bool {
//		return row >= 0 && row < rows && col >= 0 && col < columns
//	}
//	/// Returns true if input dimensions are valid.
//	public func isValidMatrix() -> Bool {
//		rows*columns == bufferSize && rows > 0 && columns > 0
//	}
////	/// Fetches or updates a value in the matrix.
////	public subscript(safely r: Int, c: Int) -> T? {
////		get {
////			guard isValidIndex(row: r, col: c) else {
////				return nil
////				print("Matrix index out of bounds: (\(r), \(c))")
////			}
////			return flatValues[unchecked: r * columns + c]
////		}
////		set(newValue) {
////			guard isValidIndex(row: r, col: c) else {
////				print("Matrix index out of bounds: (\(r), \(c))")
////			}
////			if let newValue {
////				flatValues[unchecked: r * columns + c] = newValue
////			}
////		}
////	}
//	/// Unsafely fetches or updates a value in the matrix.
//	private subscript(r: Int, c: Int) -> T {
//		get {
//			return flatValues[unchecked: r * columns + c]
//		}
//		set(newValue) {
//			flatValues[unchecked: r * columns + c] = newValue
//		}
//	}
//	
//	// MARK: Surge
//	// Surge MIT License:
//
//	//Copyright © 2014-2019 the Surge contributors
//	//
//	//Permission is hereby granted, free of charge, to any person obtaining a copy
//	//of this software and associated documentation files (the "Software"), to deal
//	//in the Software without restriction, including without limitation the rights
//	//to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
//	//copies of the Software, and to permit persons to whom the Software is
//	//furnished to do so, subject to the following conditions:
//	//
//	//The above copyright notice and this permission notice shall be included in
//	//all copies or substantial portions of the Software.
//	//
//	//THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//	//IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//	//FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
//	//AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//	//LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
//	//OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
//	//THE SOFTWARE.
//	
//	public func eigenvalues(maxIterations: Int = 100, tolerance: Double = 1e-10) throws -> [(x: T, i: T)] {
//		if rows == 1 {
//			return self[0,0].isNaN ? [] : [((self[0,0]), T(0))]
//		}
//		if rows == 2 {
//			let a = self[0,0]
//			let b = self[0,1]
//			let c = self[1,0]
//			let d = self[1,1]
//
//			let trace = a + d
//			let determinant = a * d - b * c
//			let discriminant = trace * trace - 4 * determinant
//
//			if discriminant >= 0 {
//				let sqrtDiscriminant = discriminant.squareRoot()
//				let lambda1 = (trace + sqrtDiscriminant) / 2
//				let lambda2 = (trace - sqrtDiscriminant) / 2
//				return [(lambda1,0.0), (lambda2,0.0)]
//			}
//		}
//		return try eigenDecompose(computeEigenVectors: false).eigenValues
//	}
//	
//	/// Decomposes a square matrix into its eigenvalues and left and right eigenvectors.
//	/// The decomposition may result in complex numbers, represented by (T, T), which
//	///   are the (real, imaginary) parts of the complex number.
//	/// - Returns: a struct with the eigen values and left and right eigen vectors using (T, T)
//	///   to represent a complex number.
//	func eigenDecompose(computeEigenVectors: Bool) throws -> MatrixEigenDecompositionResult<T> {
//		let input = InlineSquareMatrix<n, bufferSize, Double>(self)
//		let decomposition = try eigenDecomposeOLD(input, computeEigenVectors: computeEigenVectors)
//		
////		guard rows == columns else {
////			throw MatrixError.notSquareMatrix
////		}
////		let (Q,R) = self.qrDecompose()
////		let lambdas: [(x: T, i: T)] = []
////		var rightEigenVectors: [[(T, T)]] = [[]]
////		if computeEigenVectors  {
////			for lambda in lambdas {
////				var A = try self-Matrix<T>.identity(size: rows)*lambda.x
////				let realVector = try Matrix.solve(A: A, b: Matrix(Array<T>(repeating: 0, count: rows), isRow: false))
////				A = try self-Matrix<T>.identity(size: rows)*lambda.i
////				let imaginaryVector = try Matrix.solve(A: A, b: Matrix(Array<T>(repeating: 0, count: rows), isRow: false))
////				var eigenVector: [(T, T)] = []
////				for i in 0..<rows {
////					eigenVector.append((realVector.flatValues[i],imaginaryVector.flatValues[i]))
////				}
////				rightEigenVectors.append(eigenVector)
////			}
////		}
//		
//		return MatrixEigenDecompositionResult<T>(
//			eigenValues: decomposition.eigenValues.map { (T($0.0), T($0.1)) },
//			leftEigenVectors: decomposition.leftEigenVectors.map { $0.map { (T($0.0), T($0.1)) } },
//			rightEigenVectors: decomposition.rightEigenVectors.map { $0.map { (T($0.0), T($0.1)) } }
//		)
//	}
//
//	/// Decomposes a square matrix into its eigenvalues and left and right eigenvectors.
//	/// The decomposition may result in complex numbers, represented by (Double, Double), which
//	///   are the (real, imaginary) parts of the complex number.
//	/// - Parameters:
//	///   - lhs: a square matrix
//	/// - Returns: a struct with the eigen values and left and right eigen vectors using (Double, Double)
//	///   to represent a complex number.
//	func eigenDecomposeOLD(_ lhs: InlineSquareMatrix<n, bufferSize, Double>, computeEigenVectors: Bool) throws -> MatrixEigenDecompositionResult<Double> {
//
//		// dgeev_ needs column-major matrices, so transpose 'lhs'.
//		var matrixGrid: [__CLPK_doublereal] = Array(lhs.transpose().flatValues)
//		var matrixRowCount = __CLPK_integer(lhs.rows)
//		let matrixColCount = matrixRowCount
//		var eigenValueCount = matrixRowCount
//		var leftEigenVectorCount = matrixRowCount
//		var rightEigenVectorCount = matrixRowCount
//
//		var workspaceQuery: Double = 0.0
//		var workspaceSize = __CLPK_integer(-1)
//		var error: __CLPK_integer = 0
//
//		var eigenValueRealParts = [Double](repeating: 0, count: Int(eigenValueCount))
//		var eigenValueImaginaryParts = [Double](repeating: 0, count: Int(eigenValueCount))
//		var leftEigenVectorWork = [Double](repeating: 0, count: Int(leftEigenVectorCount * matrixColCount))
//		var rightEigenVectorWork = [Double](repeating: 0, count: Int(rightEigenVectorCount * matrixColCount))
//
//		var decompositionJobVL: [CChar]
//		var decompositionJobVR: [CChar]
//		
//		if computeEigenVectors  {
//			decompositionJobVL = [0x56, 0x00] // "V" (compute)
//			decompositionJobVR = [0x56, 0x00] // "V" (compute)
//		} else {
//			decompositionJobVL = Array("N".utf8CString) // "N" (do not compute)
//			decompositionJobVR = Array("N".utf8CString) // "N" (do not compute)
//		}
//
//		// Call dgeev to find out how much workspace to allocate
//		dgeev_(&decompositionJobVL, &decompositionJobVR, &matrixRowCount, &matrixGrid, &eigenValueCount, &eigenValueRealParts, &eigenValueImaginaryParts, &leftEigenVectorWork, &leftEigenVectorCount, &rightEigenVectorWork, &rightEigenVectorCount, &workspaceQuery, &workspaceSize, &error)
//		if error != 0 {
//			throw MatrixError.matrixNotDecomposable
//		}
//
//		// Allocate the workspace and call dgeev again to do the actual decomposition
//		var workspace = [Double](repeating: 0.0, count: Int(workspaceQuery))
//		workspaceSize = __CLPK_integer(workspaceQuery)
//		dgeev_(&decompositionJobVL, &decompositionJobVR, &matrixRowCount, &matrixGrid, &eigenValueCount, &eigenValueRealParts, &eigenValueImaginaryParts, &leftEigenVectorWork, &leftEigenVectorCount, &rightEigenVectorWork, &rightEigenVectorCount, &workspace, &workspaceSize, &error)
//		if error != 0 {
//			throw MatrixError.matrixNotDecomposable
//		}
//
//		return MatrixEigenDecompositionResult<Double>(rowCount: lhs.rows, eigenValueRealParts: eigenValueRealParts, eigenValueImaginaryParts: eigenValueImaginaryParts, leftEigenVectorWork: leftEigenVectorWork, rightEigenVectorWork: rightEigenVectorWork)
//	}
//	
//	/// Holds the result of eigendecomposition. The (Scalar, Scalar) used
//	/// in the property types represents a complex number with (real, imaginary) parts.
//	struct MatrixEigenDecompositionResult<Scalar: BinaryFloatingPoint> {
//		public let eigenValues: [(Scalar, Scalar)]
//		public let leftEigenVectors: [[(Scalar, Scalar)]]
//		public let rightEigenVectors: [[(Scalar, Scalar)]]
//		
//		public init(eigenValues: [(Scalar, Scalar)], leftEigenVectors: [[(Scalar, Scalar)]], rightEigenVectors: [[(Scalar, Scalar)]]) {
//			self.eigenValues = eigenValues
//			self.leftEigenVectors = leftEigenVectors
//			self.rightEigenVectors = rightEigenVectors
//		}
//		
//		public init(rowCount: Int, eigenValueRealParts: [Scalar], eigenValueImaginaryParts: [Scalar], leftEigenVectorWork: [Scalar], rightEigenVectorWork: [Scalar]) {
//			// The eigenvalues are an array of (real, imaginary) results from dgeev
//			self.eigenValues = Array(zip(eigenValueRealParts, eigenValueImaginaryParts))
//			
//			// Build the left and right eigenvectors
//			let emptyVector = [(Scalar, Scalar)](repeating: (0.0, 0.0), count: rowCount)
//			var leftEigenVectors = [[(Scalar, Scalar)]](repeating: emptyVector, count: rowCount)
//			MatrixEigenDecompositionResult.buildEigenVector(eigenValueImaginaryParts: eigenValueImaginaryParts, eigenVectorWork: leftEigenVectorWork, result: &leftEigenVectors)
//			
//			var rightEigenVectors = [[(Scalar, Scalar)]](repeating: emptyVector, count: rowCount)
//			MatrixEigenDecompositionResult.buildEigenVector(eigenValueImaginaryParts: eigenValueImaginaryParts, eigenVectorWork: rightEigenVectorWork, result: &rightEigenVectors)
//			
//			self.leftEigenVectors = leftEigenVectors
//			self.rightEigenVectors = rightEigenVectors
//		}
//		
//		// Convert the result of dgeev into an array of complex numbers
//		// See Intel's documentation on column-major results for sample code that this is based on:
//		// https://software.intel.com/sites/products/documentation/doclib/mkl_sa/11/mkl_lapack_examples/dgeev.htm
//		private static func buildEigenVector<S>(eigenValueImaginaryParts: [S], eigenVectorWork: [S], result: inout [[(S, S)]]) where S: FloatingPoint & ExpressibleByFloatLiteral {
//			// row and col count are the same because result must be square.
//			let rowColCount = result.count
//
//			for row in 0..<rowColCount {
//				var col = 0
//				while col < rowColCount {
//					if eigenValueImaginaryParts[col] == 0.0 {
//						// v is column-major
//						result[row][col] = (eigenVectorWork[row + rowColCount * col], 0.0)
//						col += 1
//					} else {
//						// v is column-major
//						result[row][col] = (eigenVectorWork[row + col * rowColCount], eigenVectorWork[row + rowColCount * (col + 1)])
//						result[row][col + 1] = (eigenVectorWork[row + col * rowColCount], -eigenVectorWork[row + rowColCount * (col + 1)])
//						col += 2
//					}
//				}
//			}
//		}
//
//	}
//}


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



