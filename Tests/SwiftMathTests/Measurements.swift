//
//  MeasurementTests.swift
//
//
//  Created by Jacob W Esselstyn on 6/23/23.
//

import XCTest
import MySwift
import SwiftMath

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
}
