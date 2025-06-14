//
//  MeasurementTests.swift
//
//
//  Created by Jacob W Esselstyn on 6/23/23.
//

import XCTest
import MySwift
import MySwiftMath

final class MeasurementTests: XCTestCase {
	
//	override func setUpWithError() throws {
//		// Put setup code here. This method is called before the invocation of each test method in the class.
//	}
//
//	override func tearDownWithError() throws {
//		// Put teardown code here. This method is called after the invocation of each test method in the class.
//	}

	func testInvertedTemperature() {
		
		// Use XCTAssert and related functions to verify your tests produce the correct results.
		// Any test you write for XCTest can be annotated as throws and async.
		// Mark your test throws to produce an unexpected failure when your test encounters an uncaught error.
		// Mark your test async to allow awaiting for asynchronous code to complete. Check the results with assertions afterwards.
		
//		self.measure {
//			// Put the code you want to measure the time of here.
//		}
		
		let a1 = Measurement<UnitInverseTemperature>(value: 61, unit: .inverseCelsius)
		let a1v = a1.converted(to: .inverseFahrenheit).value
		let a2 = Measurement<UnitInverseTemperature>(value: a1v, unit: .inverseFahrenheit).converted(to: .inverseCelsius)
		XCTAssertEqual(a1.value, a2.value, accuracy: 0.001)
		
		let b1 = Measurement<UnitInverseTemperature>(value: 0, unit: .inverseFahrenheit)
		let b1v = b1.converted(to: .inverseCelsius).value
		let b2 = Measurement<UnitInverseTemperature>(value: b1v, unit: .inverseCelsius).converted(to: .inverseFahrenheit)
		XCTAssertEqual(b1.value, b2.value, accuracy: 0.01)
		
		let c1 = Measurement<UnitInverseTemperature>(value: -40, unit: .inverseKelvin)
		let c1v = c1.converted(to: .inverseCelsius).value
		let c2 = Measurement<UnitInverseTemperature>(value: c1v, unit: .inverseCelsius).converted(to: .inverseKelvin)
		XCTAssertEqual(c1.value, c2.value, accuracy: 0.001)
		
		let d1 = Measurement<UnitInverseTemperature>(value: 273.15, unit: .inverseKelvin)
		let d1v = d1.converted(to: .inverseFahrenheit).value
		let d2 = Measurement<UnitInverseTemperature>(value: d1v, unit: .inverseFahrenheit).converted(to: .inverseKelvin)
		XCTAssertEqual(d1.value, d2.value, accuracy: 0.001)
		
		let e1 = Measurement<UnitInverseTemperature>(value: 32, unit: .inverseFahrenheit)
		let e1v = e1.converted(to: .inverseKelvin).value
		let e2 = Measurement<UnitInverseTemperature>(value: e1v, unit: .inverseKelvin).converted(to: .inverseFahrenheit)
		XCTAssertEqual(e1.value, e2.value, accuracy: 0.001)
	}
	
	func testUnitPressureConversions() {
		var a1 = Measurement<UnitPressure>(value: 29, unit: .kipsPerSquareInch)
		a1.convert(to: .megapascals)
		let a2 = Measurement<UnitPressure>(value: 199.9479615, unit: .megapascals)
		XCTAssertEqual(a1.value, a2.value, accuracy: 0.001)
		XCTAssertEqual(a1.unit, a2.unit)
	}
	func testUnitWorkConversions() {
		var a1 = Measurement<UnitWork>(value: 1, unit: .newtonMeters)
		a1.convert(to: .poundFeet)
		XCTAssertEqual(a1.value, 0.7375610332, accuracy: 0.001)
	}
	
	func testUnitAngleConversions() {
		var a1 = Measurement<UnitAngle>(value: 45, unit: .degrees)
		a1.convert(to: .slope)
		XCTAssertEqual(a1.value, 1, accuracy: 0.00000000001)
		
		a1 = Measurement<UnitAngle>(value: 100, unit: .percentSlope)
		a1.convert(to: .degrees)
		XCTAssertEqual(a1.value, 45, accuracy: 0.00000000001)
		
		a1 = Measurement<UnitAngle>(value: Double.infinity, unit: .percentSlope)
		a1.convert(to: .degrees)
		XCTAssertEqual(a1.value, 90, accuracy: 0.00000000001)
		
		a1 = Measurement<UnitAngle>(value: -Double.infinity, unit: .slope)
		a1.convert(to: .degrees)
		XCTAssertEqual(a1.value, -90, accuracy: 0.00000000001)
		
		a1 = Measurement<UnitAngle>(value: 90, unit: .degrees)
		a1.convert(to: .slope)
		XCTAssertEqual(a1.value, Double.infinity, accuracy: 0.00000000001)
		
		a1 = Measurement<UnitAngle>(value: -90, unit: .degrees)
		a1.convert(to: .percentSlope)
		XCTAssertEqual(a1.value, -Double.infinity, accuracy: 0.00000000001)
	}
	
	func testUnitForceConversions() {
		var a1 = Measurement<UnitForce>(value: 1, unit: .kip)
		a1.convert(to: .pound)
		XCTAssertEqual(a1.value, 1000, accuracy: 0.00000001)
		a1.convert(to: .newton)
		XCTAssertEqual(a1.value, 4448.22162825, accuracy: 0.0001)
		a1.convert(to: .kilonewton)
		XCTAssertEqual(a1.value, 4.44822162825, accuracy: 0.0001)
		
		var x1 = Measurement<UnitLength>(value: 1, unit: .feet)
		x1.convert(to: .meters)
		var work = a1*x1
		work.convert(to: .kilonewtonMeters)
		work.convert(to: .kipInches)
		XCTAssertEqual(work.value, 12, accuracy: 0.0000001)
		var distributedForce = Measurement<UnitLinearForce>(value: 1, unit: .kipsPerFoot)
		distributedForce.convert(to: .kilonewtonsPerMillimeter)
		var a2 = distributedForce*x1
		a2.convert(to: .kip)
		XCTAssertEqual(a2.value, 1, accuracy: 0.00000001)
	}
}
