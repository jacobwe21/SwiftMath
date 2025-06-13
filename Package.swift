// swift-tools-version: 6.2
// The swift-tools-version declares the minimum version of Swift required to build this package.

import PackageDescription

let package = Package(
    name: "MySwiftMath",
	defaultLocalization: "en",
	platforms: [.iOS(.v18),.macOS(.v15),.macCatalyst(.v18)],
    products: [
        // Products define the executables and libraries a package produces, making them visible to other packages.
        .library(
            name: "MySwiftMath",
            targets: ["MySwiftMath"]),
    ],
	dependencies: [
		// Dependencies declare other packages that this package depends on.
		// .package(url: /* package url */, from: "1.0.0"),
		.package(url: "https://github.com/jacobwe21/MySwift", branch: "main"),
		//.package(url: "https://github.com/Jounce/Surge.git", .upToNextMajor(from: "2.3.2")),
	],
    targets: [
        // Targets are the basic building blocks of a package, defining a module or a test suite.
        // Targets can depend on other targets in this package and products from dependencies.
		.target(
			name: "MySwiftMath",
			dependencies: ["MySwift"],
			swiftSettings: [
				.define("ACCELERATE_NEW_LAPACK")  // ⬅️ Important
			]
		),
		.testTarget(
            name: "MySwiftMathTests",
            dependencies: ["MySwiftMath","MySwift"]),
    ]
		
)
