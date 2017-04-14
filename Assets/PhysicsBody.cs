using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class PhysicsBody : MonoBehaviour {
    //Center Of Mass
    Vector2? _centerOfMass;
    public Vector2 CenterOfMass {
        get {
            if (!_centerOfMass.HasValue) {
                _centerOfMass = calculateCenterOfMass();
            }
            return _centerOfMass.Value;
        }
    }
    
    //Mass
    float? _mass;
    public float Mass {
        get {
            if (!_mass.HasValue) {
                _mass = calculateMass();
            }
            return _mass.Value;
        }
    }

    //Moment Of Inertia
    float? _momentOfInertia;
    public float MomentOfInertia {
        get {
            if (!_momentOfInertia.HasValue) {
                _momentOfInertia = calculateMomentOfInertia();
            }
            return _momentOfInertia.Value;
        }
    }

    public float CurrentAngularVelocity;
    public float AngularDisplacement;
    public float CurrentAngularAcceleration;
    public Vector2 CurrentVelocity;

	void Start () {
		print(string.Format("\nMass: {0}", Mass));
		print(string.Format("\nCenterOfMass: {0}", CenterOfMass));
		print(string.Format("\nMomentOfInertia: {0}", MomentOfInertia));

        CurrentVelocity = Vector2.zero;
        CurrentAngularVelocity = 0f;
        AngularDisplacement = 0f;
        CurrentAngularAcceleration = 0f;
	}

    void FixedUpdate() {
        transform.position += (Vector3)CurrentVelocity * Time.deltaTime;
        var degreesToRotate = -CurrentAngularVelocity * Time.deltaTime * Mathf.Rad2Deg;
        //print(degreesToRotate);
		transform.RotateAround(CenterOfMass, Vector3.forward, degreesToRotate);
		AngularDisplacement += CurrentAngularVelocity * Time.deltaTime;
        _centerOfMass = null;
    }

    //applyForce(thrusterForward.position, transform.forward, 10000);
    public void ApplyForce(Vector2 position, float force) {
        var direction = new Vector2(transform.right.y, transform.right.x).normalized;
        print( transform.eulerAngles.z + " " + direction);
        var F = force; //Newtons of force
        var R = Vector2.Distance(CenterOfMass, position); //Distance of force position and COM
        var forceDirectionWorld = position + direction;
        var comDifference = CenterOfMass - position; //Difference of force position and COM

        //? Two?
        //var T = Vector2.Angle (forceDirectionWorld, comDifference);
        var T = (CenterOfMass.x - position.x)/R;
		//var angularAcceleration = new Vector3(0.0f, R * F * T, 0.0f) / (float)MomentOfInertia;
		//var angularAcceleration = Vector3.Cross(comDifference, force * direction) / (float)MomentOfInertia;
		var angularAcceleration =  R * F * T / (float)MomentOfInertia;
		var accelerationMagnitude = F / Mass; 
		var acceleration = direction * (float)accelerationMagnitude * Time.deltaTime; // m/s

		CurrentAngularAcceleration = angularAcceleration;
		CurrentAngularVelocity = CurrentAngularVelocity + CurrentAngularAcceleration * Time.deltaTime;
		CurrentVelocity += acceleration;
    }

    #region Variable Calculation Methods

    private float calculateMass() {
        var mass = 0f;
        foreach( var child in GetComponentsInChildren<PhysicsWeight>() ) {
            mass += child.Mass;
        }
        return mass;
    }

    private Vector3 calculateCenterOfMass() {
       var com = Vector2.zero;
        foreach ( var child in GetComponentsInChildren<PhysicsWeight>() ) {
            var childCenterPoint = child.GetComponent<Renderer>().bounds.center;
            com.x += child.Mass * childCenterPoint.x;
            com.y += child.Mass * childCenterPoint.y;
        }
        com /= Mass;
        return com;
    }

    private float calculateMomentOfInertia() {
        var inertia = 0.0f;
        foreach ( var child in GetComponentsInChildren<PhysicsWeight>() ) {
            var childBounds = child.GetComponent<Renderer>().bounds;
            var boundsSize = childBounds.size;
            var childInertiaToSelf = (1.0f/12.0f)*( Mass )*(boundsSize.x*boundsSize.x + boundsSize.y*boundsSize.y + boundsSize.z*boundsSize.z);
            var dif = new Vector2( childBounds.center.x - transform.position.x, childBounds.center.y - transform.position.y);
            var distanceFromCOM = Mathf.Sqrt((dif.x * dif.x) + (dif.y * dif.y));
            var childInertia = childInertiaToSelf + child.Mass * Mathf.Pow(distanceFromCOM, 2);
            inertia += childInertia;
        }
        return inertia;
    }
    #endregion
}
