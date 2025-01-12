#include <vtkActor.h>
#include <vtkCamera.h>
#include <vtkCellPicker.h>
#include <vtkCommand.h>
#include <vtkGlyph3DMapper.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkOBJReader.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSphereSource.h>
#include <vtkSmartPointer.h>
#include <set>
#include <Eigen/Dense>
#include <vtkTextActor.h>
#include <map>
#include <vtkTextProperty.h>
#include <vtkSmartPointer.h>
#include <vtkProperty2D.h>
#include <vtkActor2D.h>
#include <vtkRenderer.h>
#include <vtkTextActor3D.h>
#include <vtkTextRenderer.h>
#include <vtkFreeTypeTools.h>

class MouseInteractorStyle : public vtkInteractorStyleTrackballCamera
{
public:
    static MouseInteractorStyle* New();
    vtkTypeMacro(MouseInteractorStyle, vtkInteractorStyleTrackballCamera);

    MouseInteractorStyle()
    {
        SelectedVertexId = -1;
        IsDragging = false;
        buttonText = nullptr;
    }

    void InitializeButton()
    {
        if (!buttonText)
        {
            buttonText = vtkSmartPointer<vtkTextActor>::New();
            buttonText->SetInput("Apply ARAP");
            buttonText->GetTextProperty()->SetColor(1.0, 1.0, 1.0);
            buttonText->GetTextProperty()->SetFontSize(16);
            buttonText->GetTextProperty()->SetJustificationToRight();
            buttonText->GetTextProperty()->SetVerticalJustificationToTop();
            
            // Use normalized viewport coordinates (0-1)
            buttonText->GetPositionCoordinate()->SetCoordinateSystemToNormalizedViewport();
            buttonText->SetPosition(0.95, 0.95);  // Place in top-right corner
            
            this->GetDefaultRenderer()->AddActor2D(buttonText);
            this->GetInteractor()->GetRenderWindow()->Render();
        }
    }

    // Add this method to be called from main after setting up the interactor
    void SetupComplete()
    {
        InitializeButton();
    }

    virtual void OnKeyPress() override
    {
        vtkRenderWindowInteractor* rwi = this->Interactor;
        std::string key = rwi->GetKeySym();
        std::string keyState = rwi->GetKeySym();  // Get the key state

        if (key == "x")
        {
            std::cout << "OnKeyPress - Key: " << key 
                      << " RepeatCount: " << rwi->GetRepeatCount()
                      << " KeyCode: " << rwi->GetKeyCode()
                      << std::endl;
            // Get current mouse position
            int* clickPos = this->GetInteractor()->GetEventPosition();

            vtkSmartPointer<vtkCellPicker> picker =
                vtkSmartPointer<vtkCellPicker>::New();
            picker->SetTolerance(0.0005);

            picker->Pick(clickPos[0], clickPos[1], 0, this->GetDefaultRenderer());
            vtkActor* actor = picker->GetActor();

            if (actor == this->MeshActor)
            {
                vtkIdType id = picker->GetPointId();
                if (id != -1)
                {
                    SelectedVertexId = id;
                    double* pos = this->MeshPolyData->GetPoint(id);
                    originalPosition[0] = pos[0];
                    originalPosition[1] = pos[1];
                    originalPosition[2] = pos[2];
                    IsDragging = true;
                    this->GetInteractor()->GetRenderWindow()->SetDesiredUpdateRate(10.0);
                    return;
                }
            }
        }
        else if (key == "Escape" && IsDragging)
        {
            std::cout << "Was dragging and abort key pressed" << std::endl;
            // Add the newPos to the modified vertices
            if (SelectedVertexId != -1)
            {
                ModifiedVertices[SelectedVertexId] = Eigen::Vector3d(newPosition[0], newPosition[1], newPosition[2]);
            }
            // Cancel vertex manipulation
            IsDragging = false;
            SelectedVertexId = -1;
        }
        
        // Don't forget to forward other keys to the parent class
        vtkInteractorStyleTrackballCamera::OnKeyPress();
    }

    virtual void OnMouseMove() override
    {
        if (IsDragging && SelectedVertexId != -1)
        {
            int* mousePos = this->GetInteractor()->GetEventPosition();
            double worldPos[4];
            this->GetDefaultRenderer()->SetDisplayPoint(mousePos[0], mousePos[1], 0);
            this->GetDefaultRenderer()->DisplayToWorld();
            this->GetDefaultRenderer()->GetWorldPoint(worldPos);
            if (worldPos[3] != 0.0)
            {
                for (int i = 0; i < 3; i++)
                {
                    worldPos[i] /= worldPos[3];
                }
                // TODO: return to this; seems to work but make sure that the math is correct
                // Get the camera's view plane normal and position
                vtkCamera* camera = this->GetDefaultRenderer()->GetActiveCamera();
                double cameraPos[3], viewPlaneNormal[3];
                camera->GetPosition(cameraPos);
                camera->GetViewPlaneNormal(viewPlaneNormal);

                // Project the mouse position onto a plane perpendicular to view
                double t = (
                    viewPlaneNormal[0] * (originalPosition[0] - worldPos[0]) +
                    viewPlaneNormal[1] * (originalPosition[1] - worldPos[1]) +
                    viewPlaneNormal[2] * (originalPosition[2] - worldPos[2])
                ) / (
                    viewPlaneNormal[0] * viewPlaneNormal[0] +
                    viewPlaneNormal[1] * viewPlaneNormal[1] + 
                    viewPlaneNormal[2] * viewPlaneNormal[2]
                );

                newPosition[0] = worldPos[0] + viewPlaneNormal[0] * t;
                newPosition[1] = worldPos[1] + viewPlaneNormal[1] * t;
                newPosition[2] = worldPos[2] + viewPlaneNormal[2] * t;

                this->MeshPolyData->GetPoints()->SetPoint(SelectedVertexId, newPosition);
                this->MeshPolyData->GetPoints()->Modified();
                this->MeshMapper->Update();
                this->GetInteractor()->GetRenderWindow()->Render();
            }
            return;
        }

        // Always call parent's OnMouseMove to maintain camera navigation
        vtkInteractorStyleTrackballCamera::OnMouseMove();
    }

    virtual void OnLeftButtonDown() override 
    {
        int* clickPos = this->GetInteractor()->GetEventPosition();
        
        // Check if ARAP button was clicked
        if (IsButtonClicked(clickPos))
        {
            if (!ModifiedVertices.empty())
            {
                ARAP();
                ModifiedVertices.clear();  // Clear modifications after ARAP
            }
            return;
        }
        
        vtkInteractorStyleTrackballCamera::OnLeftButtonDown();
    }

    void ARAP()
    {
        std::cout << "ARAP with newPos: " << newPosition[0] << ", " << newPosition[1] << ", " << newPosition[2] << std::endl;
        std::cout << "originalPos: " << originalPosition[0] << ", " << originalPosition[1] << ", " << originalPosition[2] << std::endl;

        // Compute the fitting energy by computing the l2 norm between the original position and the new position
        double Efit = std::sqrt(
            std::pow(originalPosition[0] - newPosition[0], 2) +
            std::pow(originalPosition[1] - newPosition[1], 2) +
            std::pow(originalPosition[2] - newPosition[2], 2)
        );
        std::cout << "Fitting energy: " << Efit << std::endl;

        // Unknows: 3n (rot from Ereg) + 3n (trans from Efit)
        // Constraints: 3n * valence (from Ereg) + 3c (from Efit, where c is the number of manipulated control vertices)

        // First optimize for rotations -> need to solve n SVDs, based on matrices containing the pairwise edge vectors and weights

        // Have to update the deformed edge vectors here; weights and edge vectors are already updated, either when loading the mesh or after the last ARAP step
        std::cout << "Updating deformed edge vectors" << std::endl;
        UpdateDeformedEdgeVectors();

        // Compute covariance matrices and rotation matrices for each vertex in parallel
        std::cout << "Computing covariance and rotation matrices for each vertex" << std::endl;
        std::vector<Eigen::Matrix3d> RotationMatrices(MeshPolyData->GetNumberOfPoints());
        
        #pragma omp parallel for
        for (vtkIdType i = 0; i < MeshPolyData->GetNumberOfPoints(); i++) {
            // Compute covariance matrix for vertex i
            Eigen::MatrixXd CovarianceMatrix = EdgeVectors[i] * WeightMatrices[i] * DeformedEdgeVectors[i].transpose();
            
            // Compute SVD
            Eigen::JacobiSVD<Eigen::MatrixXd> svd(CovarianceMatrix, Eigen::ComputeThinU | Eigen::ComputeThinV);
            Eigen::MatrixXd U = svd.matrixU();
            Eigen::MatrixXd V = svd.matrixV();

            // Make sure det(R) > 0 by checking if det(V * U^T) < 0
            // If so, multiply last column of U by -1
            if ((V * U.transpose()).determinant() < 0) {
                U.col(2) *= -1;
            }

            // Compute rotation matrix for vertex i
            RotationMatrices[i] = V * U.transpose();
            
        } 

        std::cout << "Rotation matrices computed" << std::endl;

        // TODO: then optimize for positions; compute system (see Eq.8 + Eq.9 in the paper) derived form the gradients,
        //  take fixed positions into account, then solve using Cholesky factorization (see TAUCS library)

        // TODO: make sure to comply with this summary: precompute coefficients wij -> prefactor system matrix (see Eq.9) -> estimate rot -> estimate trans

        // TODO: return and allow for updating the UI while computing the edge vectors and weights already, maybe on a background thread later on?
        UpdateEdgeVectors();  
        UpdateWeights();      

    }

    void SetMeshActor(vtkActor* actor)
    {
        this->MeshActor = actor;
    }

    void SetMeshPolyData(vtkPolyData* polyData)
    {
        this->MeshPolyData = polyData;
        BuildVertexNeighborhoods();
        UpdateEdgeVectors();
        UpdateWeights();
    }

    void SetMeshMapper(vtkPolyDataMapper* mapper)
    {
        this->MeshMapper = mapper;
    }

private:
    vtkActor* MeshActor;
    vtkPolyData* MeshPolyData;
    vtkPolyDataMapper* MeshMapper;
    vtkIdType SelectedVertexId;
    double originalPosition[3];
    double newPosition[3];
    bool IsDragging;
    std::vector<std::set<vtkIdType>> VertexNeighbors;
    std::vector<Eigen::MatrixXd> EdgeVectors;           // Per-vertex edge matrices (3 × valence)
    std::vector<Eigen::MatrixXd> DeformedEdgeVectors;   // Same structure as EdgeVectors
    std::vector<Eigen::DiagonalMatrix<double, Eigen::Dynamic>> WeightMatrices;  // Per-vertex weights
    std::map<vtkIdType, Eigen::Vector3d> ModifiedVertices;
    vtkSmartPointer<vtkActor> ArapButton;
    vtkSmartPointer<vtkTextActor> buttonText;
    
    void BuildVertexNeighborhoods()
    {
        if (!MeshPolyData) return;

        vtkIdType numPoints = MeshPolyData->GetNumberOfPoints();
        std::cout << "Number of vertices in mesh: " << numPoints << std::endl;
        
        VertexNeighbors.resize(numPoints);
        EdgeVectors.resize(numPoints);
        DeformedEdgeVectors.resize(numPoints);
        WeightMatrices.resize(numPoints);
        
        // Build neighborhoods
        vtkCellArray* cells = MeshPolyData->GetPolys();
        vtkIdType npts;
        vtkIdType* pts;
        cells->InitTraversal();
        while (cells->GetNextCell(npts, pts))
        {
            for (vtkIdType i = 0; i < npts; i++)
            {
                vtkIdType v1 = pts[i];
                for (vtkIdType j = 0; j < npts; j++)
                {
                    if (i != j)
                    {
                        VertexNeighbors[v1].insert(pts[j]);
                    }
                }
            }
        }

        // Initialize matrices with correct sizes
        for (vtkIdType i = 0; i < numPoints; i++)
        {
            size_t valence = VertexNeighbors[i].size();
            EdgeVectors[i].resize(3, valence);
            DeformedEdgeVectors[i].resize(3, valence);
            WeightMatrices[i].resize(valence);
        }

    }

    void UpdateEdgeVectors()
    {
        for (vtkIdType i = 0; i < MeshPolyData->GetNumberOfPoints(); i++)
        {
            const auto& neighbors = VertexNeighbors[i];
            const Eigen::Vector3d pi(MeshPolyData->GetPoint(i));
            
            int col = 0;
            for (vtkIdType neighborId : neighbors)
            {
                const Eigen::Vector3d pj(MeshPolyData->GetPoint(neighborId));
                EdgeVectors[i].col(col++) = pj - pi;
            }
        }
    }

    bool IsButtonClicked(int* pos)
    {
        // Simple button hit test
        int* size = this->GetInteractor()->GetRenderWindow()->GetSize();
        return (pos[0] > size[0] - 100 && pos[0] < size[0] - 20 &&
                pos[1] > size[1] - 30 && pos[1] < size[1] - 10);
    }

    void UpdateDeformedEdgeVectors()
    {
        for (vtkIdType i = 0; i < MeshPolyData->GetNumberOfPoints(); i++)
        {
            const auto& neighbors = VertexNeighbors[i];
            
            // Get current position (original or modified)
            Eigen::Vector3d pi = (ModifiedVertices.count(i) > 0) ? 
                ModifiedVertices[i] : 
                Eigen::Vector3d(MeshPolyData->GetPoint(i));
            
            int col = 0;
            for (vtkIdType neighborId : neighbors)
            {
                // Get neighbor position (original or modified)
                Eigen::Vector3d pj = (ModifiedVertices.count(neighborId) > 0) ? 
                    ModifiedVertices[neighborId] : 
                    Eigen::Vector3d(MeshPolyData->GetPoint(neighborId));
                
                DeformedEdgeVectors[i].col(col++) = pj - pi;
            }
        }
    }

    // TODO: return to this, want to use non-uniform weights later
    void UpdateWeights()
    {
        for (vtkIdType i = 0; i < MeshPolyData->GetNumberOfPoints(); i++)
        {
            // Set uniform weights for each vertex's edges
            WeightMatrices[i].diagonal().setConstant(1.0 / VertexNeighbors[i].size());
        }
    }

};

vtkStandardNewMacro(MouseInteractorStyle);

int main(int argc, char* argv[])
{
    if (argc < 2)
    {
        std::cerr << "Usage: " << argv[0] << " Filename.obj" << std::endl;
        return EXIT_FAILURE;
    }

    // Read the OBJ file
    vtkSmartPointer<vtkOBJReader> reader =
        vtkSmartPointer<vtkOBJReader>::New();
    reader->SetFileName(argv[1]);
    reader->Update();

    // Create mapper and actor for the mesh
    vtkSmartPointer<vtkPolyDataMapper> mapper =
        vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputConnection(reader->GetOutputPort());

    vtkSmartPointer<vtkActor> meshActor =
        vtkSmartPointer<vtkActor>::New();
    meshActor->SetMapper(mapper);
    meshActor->GetProperty()->SetColor(0.8, 0.8, 0.8);

    // Create sphere glyphs for vertices
    vtkSmartPointer<vtkSphereSource> sphereSource =
        vtkSmartPointer<vtkSphereSource>::New();
    sphereSource->SetRadius(0.001);
    sphereSource->Update();

    vtkSmartPointer<vtkGlyph3DMapper> glyphMapper =
        vtkSmartPointer<vtkGlyph3DMapper>::New();
    glyphMapper->SetInputData(reader->GetOutput());
    glyphMapper->SetSourceConnection(sphereSource->GetOutputPort());
    glyphMapper->Update();

    vtkSmartPointer<vtkActor> glyphActor =
        vtkSmartPointer<vtkActor>::New();
    glyphActor->SetMapper(glyphMapper);
    glyphActor->GetProperty()->SetColor(1.0, 0.0, 0.0);

    // Renderer setup
    vtkSmartPointer<vtkRenderer> renderer =
        vtkSmartPointer<vtkRenderer>::New();
    renderer->AddActor(meshActor);
    renderer->AddActor(glyphActor);
    renderer->SetBackground(0.1, 0.2, 0.4);

    // Render window
    vtkSmartPointer<vtkRenderWindow> renderWindow =
        vtkSmartPointer<vtkRenderWindow>::New();
    renderWindow->AddRenderer(renderer);
    renderWindow->SetWindowName("Interactive OBJ Viewer");

    // Interactor
    vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
        vtkSmartPointer<vtkRenderWindowInteractor>::New();
    renderWindowInteractor->SetRenderWindow(renderWindow);

    // Set custom interactor style
    vtkSmartPointer<MouseInteractorStyle> style =
        vtkSmartPointer<MouseInteractorStyle>::New();
    style->SetDefaultRenderer(renderer);
    style->SetMeshActor(meshActor);
    style->SetMeshPolyData(reader->GetOutput());
    style->SetMeshMapper(mapper);
    renderWindowInteractor->SetInteractorStyle(style);
    renderWindowInteractor->Initialize();
    
    // Add this line after interactor initialization
    style->SetupComplete();
    
    // Start interaction
    renderWindow->Render();
    renderWindowInteractor->Start();

    return EXIT_SUCCESS;
}