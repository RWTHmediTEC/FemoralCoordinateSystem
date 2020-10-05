function props = inertiaInfo(mesh)

% Get Volume (V), Center of Mass (CoM), Inertia Tensor (J) of the Bone
[props.V, props.CoM, props.J] = VolumeIntegrate(mesh.vertices, mesh.faces);

% Get Principal Axes (pAxes) & Principal Moments of Inertia (Jii)
% Sign of the Eigenvectors can change (in agreement with their general definition)
[props.pAxes, props.Jii] = eig(props.J);

% Keep the determinant positive
if det(props.pAxes) < 0
    props.pAxes = -1*props.pAxes;
end

% Create a transformation to move the mesh in the coordinate system based 
% on the principial axes of inertia and the center of mass of the mesh.
props.InertiaTFM = [props.pAxes' [0 0 0]'; 0 0 0 1]*[eye(3) -props.CoM; 0 0 0 1];

end