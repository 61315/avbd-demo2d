/*
* Copyright (c) 2025 Chris Giles
*
* Permission to use, copy, modify, distribute and sell this software
* and its documentation for any purpose is hereby granted without fee,
* provided that the above copyright notice appear in all copies.
* Chris Giles makes no representations about the suitability
* of this software for any purpose.
* It is provided "as is" without express or implied warranty.
*/

#include "solver.h"

Spring::Spring(Solver* solver, Rigid* bodyA, Rigid* bodyB, float2 rA, float2 rB, float stiffness, float rest)
    : Force(solver, bodyA, bodyB), rA(rA), rB(rB), rest(rest)
{
    this->stiffness[0] = stiffness;
    if (this->rest < 0)
        this->rest = length(transform(bodyA->position, rA) - transform(bodyB->position, rB));
}

void Spring::computeConstraint(float alpha)
{
    // Compute constraint function at current state C(x)
    C[0] = length(transform(bodyA->position, rA) - transform(bodyB->position, rB)) - rest;
}

void Spring::computeDerivatives(Rigid* body)
{
    // Compute the first and second derivatives for the desired body
    // R = [c -s; s c]; and dR/dtheta = S = [-s -c; c -s]; 
    // if theta = 0, R = [1 0; 0 1]; and S = [0 -1; 1 0];
    float2x2 S = { 0, -1, 1, 0 };
    float2x2 I = { 1, 0, 0, 1 };

    // .transform(), point transform
    // .rotate(), vector transform

    // attachment point delta in world coords
    float2 d = transform(bodyA->position, rA) - transform(bodyB->position, rB);
    float dlen2 = dot(d, d);
    if (dlen2 == 0)
        return;
    float dlen = sqrtf(dlen2);
    float2 n = d / dlen; // unit direction in world coords (rB_world to rA_world)
    
    // double partial of `d` wrt `x`
    float2x2 dxx = (I - outer(n, n) / dlen2) / dlen;

    // For body_a:
    // rAw = R(theta_a) * rA
    // d = rAw - rBw
    // |d| = |rAw - rBw|
    // C = |d| - rest
    // S = dR(theta_a)/dtheta_a = a 90-degree turn ccw
    // dd/dtheta_a = drAw/dtheta_a = S * R(theta_a) * rA = S * rAw
    // dx = dC/dx_a = d/|d| or `n`
    // dr = dC/dtheta_a = d|d|/dtheta_a = n * (dd/dtheta_a) = n * (S * rAw)
    // dxx = dn/dx = (I - n*n')/|d|
    // dxr = dn/dtheta_a = dxx * S * rAw
    //            ^--------------------------->
    // drr = d(n * S * rAw)/dtheta_a = (dn/dtheta_a) * S * rAw + n * S * (drAw/dtheta_a)
    //                               = (dxr) * S * rAw         + n * S * (S * rAw)
    //                               = dxr * S * rAw           + n * (S * S) * rAw
    //                               = dxr * S * rAw           + n * -I * rAw
    //                               = dot(dxr, S * rAw))      - dot(n, rAw)
    if (body == bodyA)
    {
        // theta is stored in position.z()
        float2 Sr = rotate(bodyA->position.z, S * rA); // Sr = S * R(theta_a) * rA = S * rAw
        float2 r = rotate(bodyA->position.z, rA); // r = R(theta_a) * rA, or call it rAw
        float2 dxr = dxx * Sr; // dxx * S * rAw
        float drr = dot(Sr, dxr) - dot(n, r);

        J[0].xy() = n; // dC/dx_a = dx = n
        J[0].z = dot(n, Sr); // dC/dtheta_a = dr = n * S * rAw 
        H[0] = { // assemble with dxx, dxr, drr
            dxx.row[0].x, dxx.row[0].y, dxr.x,
            dxx.row[1].x, dxx.row[1].y, dxr.y,
            dxr.x,          dxr.y,        drr
        };
    }
    else
    {
        float2 Sr = rotate(bodyB->position.z, S * rB);
        float2 r = rotate(bodyB->position.z, rB);
        float2 dxr = dxx * -Sr;
        float drr = dot(Sr, dxr) + dot(n, r);

        J[0].xy() = -n;
        J[0].z = dot(n, -Sr);
        H[0] = {
            dxx.row[0].x, dxx.row[0].y, dxr.x,
            dxx.row[1].x, dxx.row[1].y, dxr.y,
            dxr.x,          dxr.y,        drr
        };
    }
}

void Spring::draw() const
{
    float2 v0 = transform(bodyA->position, rA);
    float2 v1 = transform(bodyB->position, rB);

    glColor3f(0.75f, 0.0f, 0.0f);
    glBegin(GL_LINES);
    glVertex2f(v0.x, v0.y);
    glVertex2f(v1.x, v1.y);
    glEnd();
}
