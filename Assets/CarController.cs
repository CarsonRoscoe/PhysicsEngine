using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class CarController : MonoBehaviour {
    private PhysicsBody physicsBody;
    private Transform thrusterForward;
    private Transform thrusterTurnRight;
    private Transform thrusterTurnLeft;

	
	void Start () {
		physicsBody = GetComponent<PhysicsBody>();
        thrusterForward = GameObject.Find("ThrusterPoint").transform;
        thrusterTurnRight = GameObject.Find("RightTurnPoint").transform;
        thrusterTurnLeft = GameObject.Find("LeftTurnPoint").transform;
	}
	
	void Update () {
        //Gas
		if (Input.GetKey(KeyCode.W)) {

            physicsBody.ApplyForce(thrusterForward.position, 15000);
        }

        //Break/reverse
        if (Input.GetKey(KeyCode.S)) {
            physicsBody.ApplyForce(thrusterForward.position, -15000);
        }

        //Turn left
        if (Input.GetKey(KeyCode.A)) {
            physicsBody.ApplyForce(thrusterTurnLeft.position, 10000);
        }

        //Turn right
        if (Input.GetKey(KeyCode.D)) {
            physicsBody.ApplyForce(thrusterTurnRight.position, 10000);
        }
	}
}
