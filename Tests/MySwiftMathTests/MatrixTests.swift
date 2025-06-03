//
//  MatrixTests.swift
//  MySwiftMath
//
//  Created by Jacob W Esselstyn on 6/1/25.
//

import XCTest
import MySwift
import MySwiftMath

final class MatrixTests: XCTestCase {
	
//	override func setUpWithError() throws {
//		// Put setup code here. This method is called before the invocation of each test method in the class.
//	}
//
//	override func tearDownWithError() throws {
//		// Put teardown code here. This method is called after the invocation of each test method in the class.
//	}

	func testPermutationMatrix() {
		let P = [2, 0, 1]
		let Pmat = Matrix<Double>.permutationMatrix(from: P)

		for row in Pmat.values {
			print(row)
		}
		let Pmatrix = Matrix<Double>([
			[0,1,0],
			[0,0,1],
			[1,0,0]
		])
		XCTAssertTrue(Pmatrix ~= Pmat)
	}
	
	func testDeterminant() {
		var A = Matrix<Double>([
			[2,3],
			[4,2],
		])
		var det = try! A.determinant()
		XCTAssertEqual(det,-8.0)
		
		A = Matrix<Double>([
			[2,3,8],
			[4,2,3],
			[4,2,1],
		])
		det = try! A.determinant()
		XCTAssertEqual(det,16.0)
		
		A = Matrix<Double>([
			[2,3,8,5],
			[4,2,3,9],
			[4,2,1,5],
			[4,-3,1,5],
		])
		det = try! A.determinant()
		XCTAssertEqual(det,-500)
		
		A = Matrix<Double>([
			[2,3,8,5,1],
			[4,2,3,9,9],
			[4,2,1,5,2],
			[4,-3,1,5,9],
			[4,-3,1,3,-2],
		])
		det = try! A.determinant()
		XCTAssertTrue((det-3624).isApproxEqual(to: 0))
		
		A = Matrix<Double>([
			[2,3,8,5,1],
			[4,2,3,9,9],
			[4,2,1,5,2],
			[4,-3,1,5,9],
			[4,-3,1,5,9],
		])
		det = try! A.determinant()
		XCTAssertEqual(det,0)
		
		
		A = Matrix<Double>([
			[2,3,8,5,1,10,2],
			[4,2,3,9,9,19,8],
			[4,2,1,5,2,10,2],
			[4,-3,1,5,9,19,1],
			[4,-3,-4,5,2,1,9],
			[2,-3,3,8,18,3,4],
			[1,3,2,-8,13,3,2],
		])
		det = try! A.determinant()
		XCTAssertTrue((det+4715536).isApproxEqual(to: 0))
		
		A = Matrix<Double>([
			[2,3,8,5,1,10,2],
			[4,2,3,-9,9,19,8],
			[4,2,1,-4,2,10,2],
			[4,-3,1,5,9,19,1],
			[4,-3,-4,5,2,1,9],
			[2,-3,3,8,18,3,4],
			[1,3,2,-8,13,3,2],
		])
		det = try! A.determinant()
		XCTAssertTrue((det-1681988).isApproxEqual(to: 0))
	}
	
	func testLUDecomposition() {
		let n = 4
		for _ in 0...20 {
			var A = Matrix<Double>(size: n)
			for r in 0..<A.rows {
				for c in 0..<A.columns {
					A[r,c] = Double.random(in: -12...12)
				}
			}
			do {
				let (L,U,P,_) = try Matrix.luDecomposition(of: A.values)
				let pMatrix = Matrix<Double>.permutationMatrix(from: P)
				let result = try pMatrix*Matrix(L)*Matrix(U)
				print(pMatrix)
				print(result)
				print(A)
				XCTAssertTrue(result ~= A)
			} catch {
				XCTFail(error.localizedDescription)
			}
		}
		for _ in 0...20 {
			var A = Matrix<Double>(size: n)
			for r in 0..<A.rows {
				for c in 0..<A.columns {
					A[r,c] = Double.random(in: -12...12)
				}
			}
			for i in 1..<n {
				if Bool.random() {
					A[i,i] = 0.0
				}
			}
			do {
				let (L,U,P,_) = try Matrix.luDecomposition(of: A.values)
				let pMatrix = Matrix<Double>.permutationMatrix(from: P)
				let result = try pMatrix*Matrix(L)*Matrix(U)
				print(pMatrix)
				print(result)
				print(A)
				XCTAssertTrue(result ~= A)
			} catch {
				print("This matrix is singular")
			}
		}
	}
	
	func testInversion() {
		for _ in 0...20 {
			for n in 1...8 {
				var A = Matrix<Double>(size: n)
				for r in 0..<A.rows {
					for c in 0..<A.columns {
						A[r,c] = Double.random(in: 1...12)
					}
				}
				for i in 1..<n {
					if Bool.random() {
						A[i,i] = 0.0
					}
				}
				do {
					let inverse = try A.inverse()
					let identity = Matrix<Double>.identity(size: n)
					let calculatedIdentity = try A*inverse
					XCTAssertTrue(calculatedIdentity ~= identity)
				} catch {
					if let matrixError = error as? Matrix<Double>.MatrixError {
						let determinant = try? A.determinant()
						XCTAssertTrue(matrixError == .singularMatrix)
						XCTAssertTrue(determinant == 0)
						print("This matrix is singular")
					} else {
						XCTFail(error.localizedDescription)
					}
				}
			}
		}
	}
	
	func testSolve() {
		for _ in 0...20 {
			for n in 2...8 {
				var A = Matrix<Double>(size: n)
				for r in 0..<A.rows {
					for c in 0..<A.columns {
						A[r,c] = Double.random(in: -12...12)
					}
				}
				for i in 1..<n {
					if Bool.random() {
						A[i,i] = 0.0
					}
				}
				var b = Matrix<Double>(rows: n, columns: 1)
				for r in 0..<b.rows {
					b[r,0] = Double.random(in: -12...12)
				}
				do {
					let x = try Matrix.solve(A: A, b: b)
					let calculatedB = try A*x
					XCTAssertTrue(calculatedB ~= b)
				} catch {
					if let matrixError = error as? Matrix<Double>.MatrixError {
						let determinant = try? A.determinant()
						XCTAssertTrue(matrixError == .singularMatrix)
						XCTAssertTrue(determinant == 0)
						print("This matrix is singular")
					} else {
						XCTFail(error.localizedDescription)
					}
				}
			}
		}
	}
	
}
