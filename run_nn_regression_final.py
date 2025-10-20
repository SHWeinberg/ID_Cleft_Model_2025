import os
import time
import torch
import torch.nn as nn
import torch.optim as optim
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.model_selection import KFold
from sklearn.preprocessing import StandardScaler
import os
import matplotlib.pyplot as plt
import seaborn as sns


# set_seed(42)

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
print(f"Using device: {device}")

class RegressionNN(nn.Module):
    def __init__(self, input_dim, output_dim):
        super(RegressionNN, self).__init__()
        self.fc1 = nn.Linear(input_dim, 256)
        self.fc2 = nn.Linear(256, 128)
        self.fc3 = nn.Linear(128, output_dim)
        self.relu = nn.ReLU()
        self.dropout = nn.Dropout(0.0)
    def forward(self, x):
        x = self.relu(self.fc1(x))
        x = self.dropout(x)
        x = self.relu(self.fc2(x))
        return self.fc3(x)

def gradient_stability(gradient_matrixs, n_folds, output_dim):
    gradient_vectors = []
    for i in range(n_folds):
        g_i = gradient_matrixs[:, i*output_dim:(i+1)*output_dim]  # shape: (input_dim, output_dim)
        g_vec = g_i.flatten()
        gradient_vectors.append(g_vec)

    gradient_vectors = np.stack(gradient_vectors, axis=0)  # shape: (n_folds, input_dim * output_dim)
    cos_sim_matrix = cosine_similarity(gradient_vectors)
    return cos_sim_matrix


def visualize_feature_influence(gradient_matrix, input_names, output_names, title, save_dir):
    grad_matrix_np = gradient_matrix.cpu().numpy()
    abs_max = np.max(np.abs(grad_matrix_np))

    plt.figure(figsize=(18, 10))

    sns.set_theme(font_scale=1.2)

    ax = sns.heatmap(
        grad_matrix_np, 
        cmap="coolwarm", 
        annot=True, 
        fmt=".2e",
        xticklabels=output_names, 
        yticklabels=input_names,
        vmin=-abs_max, 
        vmax=abs_max, 
        center=0,
        cbar_kws={"shrink": 0.6}  
    )

    plt.title(title, fontsize=18)
    plt.xlabel("Output Features", fontsize=14)
    plt.ylabel("Input Features", fontsize=14)

    plt.tight_layout()

    filename = os.path.join(save_dir, "Influence_mat.png")
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.show()
    plt.close()
    
    print(f"Saved heatmap to: {filename}")


def analyze_feature_influence(model, sample_input):
    """
    Analyzes the influence of each input feature on the model's output
    using gradient-based analysis.

    Parameters:
    - model: The neural network model for which feature influence is analyzed.
    - sample_input: A sample input tensor to compute gradients with respect to.

    Returns:
    - gradient_matrix: A tensor containing the average gradient values
                       for each input feature with respect to each output feature.
    """
    sample_input_copy = sample_input.detach().clone()
    sample_input_copy.requires_grad = True  # Enable gradient tracking for the input
    output = model(sample_input_copy)       # Forward pass to compute model output
    num_inputs = sample_input_copy.size(1)  # Number of input features
    num_outputs = output.size(1)       # Number of output features
    gradient_matrix = torch.zeros(num_inputs, num_outputs)  # Initialize gradient matrix

    # Compute gradient for each output feature with respect to all input features
    for i in range(num_outputs):
        model.zero_grad()             # Zero out gradients before each iteration
        sample_input_copy.grad = None      # Clear existing gradients to avoid accumulation
        # Compute the gradient of the mean of the i-th output feature
        output[:, i].mean().backward(retain_graph=True)
        # Store the mean gradient for each input feature
        gradient_matrix[:, i] = sample_input_copy.grad.mean(dim=0)

    return gradient_matrix
# Loading data

start = time.time()

model_type = "Nano_CV"
D = 1
test_D_flag = 0

if model_type == "cleft":
    file_path1 = r'D:\codes\matlab code\Save_data\All_mesh_data\New_Tissue_hetg_Cleft_data_cycle1000_beats10_D1_gj_loc1_chan_loc1.csv'
    file_path2 = r'D:\codes\matlab code\Save_data\All_mesh_data\New_Tissue_hetg_Cleft_data_cycle1000_beats10_D01_gj_loc1_chan_loc1.csv'
    x_range = list(range(0, 12)) + [-1] if test_D_flag == 1 else range(0, 12)
    y_range = range(12, 24)
elif model_type == "CV":
    file_path1 = r'D:\codes\matlab code\Save_data\All_mesh_data\Tissue_hetg_property_cv_data_cycle1000_beats10_D1_gj_loc1_chan_loc1.csv'
    file_path2 = r'D:\codes\matlab code\Save_data\All_mesh_data\Tissue_hetg_property_cv_data_cycle1000_beats10_D01_gj_loc1_chan_loc1.csv'
    x_range = list(range(0, 5)) + [-1] if test_D_flag == 1 else range(0, 5)
    y_range = -2
elif model_type == "Nano_CV":
    file_path1 = r'D:\codes\matlab code\Save_data\All_mesh_data\Nano_CV_data_cycle1000_beats10_D1_gj_loc1_chan_loc1.csv'
    file_path2 = r'D:\codes\matlab code\Save_data\All_mesh_data\Nano_CV_data_cycle1000_beats10_D01_gj_loc1_chan_loc1.csv'
    x_range = list(range(0, 28)) + [-1] if test_D_flag == 1 else list(range(0, 28))
    y_range = 33
else:
    raise ValueError(f"Unknown model type: {model_type}")

df = pd.read_csv(file_path1)
df2 = pd.read_csv(file_path2)


if model_type == "cleft":
    input_name_list = df.columns[x_range].tolist()
    output_name_list = df.columns[y_range].tolist()
    if test_D_flag == 1:
        X_all = np.vstack([df.iloc[:, x_range].values, df2.iloc[:, x_range].values])
        Y_all = np.vstack([df.iloc[:, y_range].values, df2.iloc[:, y_range].values])

        base_dir = "Regression nn_multi_layer"
        main_folder = "Save_training_parameters_plots"
        sub_folder = "Cleft_All_data"
        plot_subfolder = "Plots"

        save_train_path = os.path.join(base_dir, main_folder, sub_folder)      
        save_plot_path = os.path.join(save_train_path, plot_subfolder)         
        os.makedirs(save_train_path, exist_ok=True)
        os.makedirs(save_plot_path, exist_ok=True)

    elif D == 1 and test_D_flag != 1:
        X_all = df.iloc[:, x_range].values
        Y_all = df.iloc[:, y_range].values

        base_dir = "Regression nn_multi_layer"
        main_folder = "Save_training_parameters_plots"
        sub_folder = "Cleft_D10_data"
        plot_subfolder = "Plots"

        save_train_path = os.path.join(base_dir, main_folder, sub_folder)      
        save_plot_path = os.path.join(save_train_path, plot_subfolder)         
        os.makedirs(save_train_path, exist_ok=True)
        os.makedirs(save_plot_path, exist_ok=True)
    elif D == 0.1 and test_D_flag != 1:
        base_dir = "Regression nn_multi_layer"
        main_folder = "Save_training_parameters_plots"
        sub_folder = "Cleft_D01_data"
        plot_subfolder = "Plots"

        save_train_path = os.path.join(base_dir, main_folder, sub_folder)      
        save_plot_path = os.path.join(save_train_path, plot_subfolder)         
        os.makedirs(save_train_path, exist_ok=True)
        os.makedirs(save_plot_path, exist_ok=True)
        X_all = df2.iloc[:, x_range].values
        Y_all = df2.iloc[:, y_range].values
elif model_type == "CV":
    input_name_list = df.columns[x_range].tolist()
    output_name_list = [df.columns[y_range]]
    if test_D_flag == 1:
        X_all = np.vstack([df.iloc[:, x_range].values, df2.iloc[:, x_range].values])
        Y_all = np.vstack([df.iloc[:, y_range].values.reshape(-1,1), df2.iloc[:, y_range].values.reshape(-1,1)])
        base_dir = "Regression nn_multi_layer"
        main_folder = "Save_training_parameters_plots"
        sub_folder = "CV_All_data"
        plot_subfolder = "Plots"

        save_train_path = os.path.join(base_dir, main_folder, sub_folder)      
        save_plot_path = os.path.join(save_train_path, plot_subfolder)         
        os.makedirs(save_train_path, exist_ok=True)
        os.makedirs(save_plot_path, exist_ok=True)
    elif D == 1 and test_D_flag != 1:
        X_all = df.iloc[:, x_range].values
        Y_all = df.iloc[:, y_range].values.reshape(-1,1)
        base_dir = "Regression nn_multi_layer"
        main_folder = "Save_training_parameters_plots"
        sub_folder = "CV_D10_data"
        plot_subfolder = "Plots"

        save_train_path = os.path.join(base_dir, main_folder, sub_folder)      
        save_plot_path = os.path.join(save_train_path, plot_subfolder)         
        os.makedirs(save_train_path, exist_ok=True)
        os.makedirs(save_plot_path, exist_ok=True)
    elif D == 0.1 and test_D_flag != 1:
        X_all = df2.iloc[:, x_range].values
        Y_all = df2.iloc[:, y_range].values.reshape(-1,1)

        base_dir = "Regression nn_multi_layer"
        main_folder = "Save_training_parameters_plots"
        sub_folder = "CV_D01_data"
        plot_subfolder = "Plots"

        save_train_path = os.path.join(base_dir, main_folder, sub_folder)      
        save_plot_path = os.path.join(save_train_path, plot_subfolder)         
        os.makedirs(save_train_path, exist_ok=True)
        os.makedirs(save_plot_path, exist_ok=True)
else:
    input_name_list = df.columns[x_range].tolist()
    output_name_list = [df.columns[y_range]]
    if test_D_flag == 1:
        X_all = np.vstack([df.iloc[:, x_range].values, df2.iloc[:, x_range].values])
        Y_all = np.vstack([df.iloc[:, y_range].values.reshape(-1,1), df2.iloc[:, y_range].values.reshape(-1,1)])
        base_dir = "Regression nn_multi_layer"
        main_folder = "Save_training_parameters_plots"
        sub_folder = "Nano_CV_All_data"
        plot_subfolder = "Plots"

        save_train_path = os.path.join(base_dir, main_folder, sub_folder)      
        save_plot_path = os.path.join(save_train_path, plot_subfolder)         
        os.makedirs(save_train_path, exist_ok=True)
        os.makedirs(save_plot_path, exist_ok=True)
    elif D == 1 and test_D_flag != 1:
        X_all = df.iloc[:, x_range].values
        Y_all = df.iloc[:, y_range].values.reshape(-1,1)
        base_dir = "Regression nn_multi_layer"
        main_folder = "Save_training_parameters_plots"
        sub_folder = "Nano_CV_D10_data"
        plot_subfolder = "Plots"

        save_train_path = os.path.join(base_dir, main_folder, sub_folder)      
        save_plot_path = os.path.join(save_train_path, plot_subfolder)         
        os.makedirs(save_train_path, exist_ok=True)
        os.makedirs(save_plot_path, exist_ok=True)
    elif D == 0.1 and test_D_flag != 1:
        X_all = df2.iloc[:, x_range].values
        Y_all = df2.iloc[:, y_range].values.reshape(-1,1)

        base_dir = "Regression nn_multi_layer"
        main_folder = "Save_training_parameters_plots"
        sub_folder = "Nano_CV_D01_data"
        plot_subfolder = "Plots"

        save_train_path = os.path.join(base_dir, main_folder, sub_folder)      
        save_plot_path = os.path.join(save_train_path, plot_subfolder)         
        os.makedirs(save_train_path, exist_ok=True)
        os.makedirs(save_plot_path, exist_ok=True)

scaler_X = StandardScaler()
scaler_Y = StandardScaler()
X_all = scaler_X.fit_transform(X_all)
Y_all = scaler_Y.fit_transform(Y_all)



X_tensor = torch.tensor(X_all, dtype=torch.float32).to(device)
Y_tensor = torch.tensor(Y_all, dtype=torch.float32).to(device)

n_folds = 5
kf = KFold(n_splits=n_folds, shuffle=True, random_state=None)

input_dim = X_tensor.size(1)
output_dim = Y_tensor.size(1)

criterion = nn.MSELoss()
learning_rate = 5e-4
epochs = 35
batch_size = 16

All_gradient_mat = torch.zeros(input_dim, output_dim)

output_name = output_name_list

train_loss = 0 
test_loss = 0
Y_train_plot = []
Y_test_plot  = []
Y_train_true = []
Y_test_true  = []
gradient_matrixs = []
for fold_idx, (train_idx, test_idx) in enumerate(kf.split(X_tensor)):

    X_train = X_tensor[train_idx]
    Y_train = Y_tensor[train_idx]
    X_test = X_tensor[test_idx] 
    Y_test = Y_tensor[test_idx]
    Y_train_true.append(Y_train.cpu().numpy())
    Y_test_true.append(Y_test.cpu().numpy())
    num_batches = len(X_train) // batch_size

    print(f"\n=== Training model {fold_idx + 1}/{n_folds} ===")
    current_model_path = os.path.join(save_train_path, f"KFold_model_output_{fold_idx}.pth")
    
    model = RegressionNN(input_dim, output_dim).to(device)
    
    if os.path.exists(current_model_path):
        print(f"Load model: {current_model_path}")
        model.load_state_dict(torch.load(current_model_path))
        model.eval()
    else:
        print(f"Training new model {fold_idx + 1}...")
        optimizer = optim.Adam(model.parameters(), lr=learning_rate)
        scheduler = optim.lr_scheduler.StepLR(optimizer, step_size=5, gamma=0.5)
        
        best_val_loss = float('inf')
        patience_threshold = 15 
        patience_count = 0
        best_model_state = None
        
        for epoch in range(epochs):
            model.train()
            
            perm = torch.randperm(X_train.size(0))
            X_shuf = X_train[perm]
            Y_shuf = Y_train[perm]
            
            # Train epoch
            for i in range(num_batches):
                start_idx = i * batch_size
                end_idx = min(start_idx + batch_size, X_train.size(0))
                xb = X_shuf[start_idx:end_idx]
                yb = Y_shuf[start_idx:end_idx]
                
                pred = model(xb)
                loss = criterion(pred, yb)
                
                optimizer.zero_grad()
                loss.backward()
                optimizer.step()
            
            scheduler.step()
            
            if (epoch + 1) % 5 == 0:
                model.eval()
                with torch.no_grad():
                    val_pred = model(X_test)
                    val_loss = criterion(val_pred, Y_test).item()
                
                print(f"[{output_name}] Model {fold_idx+1}, Epoch {epoch+1}/{epochs}, Test Val Loss: {val_loss:.4f}")
                
                if val_loss < best_val_loss:
                    best_val_loss = val_loss
                    patience_count = 0
                    best_model_state = model.state_dict().copy()
                    # Save best model
                    torch.save(best_model_state, current_model_path)
                    print(f"Save best model: {current_model_path}")
                else:
                    patience_count += 1
                
                if patience_count >= patience_threshold:
                    print(f"Early stopping at epoch {epoch+1}")
                    break
        
        # Load best model
        if best_model_state is not None:
            model.load_state_dict(best_model_state)
    
    print(f"Analyzing feature influence for model {fold_idx + 1}...")
    gradient_matrix = analyze_feature_influence(model, X_train)
    gradient_matrixs.append(gradient_matrix.detach().cpu().numpy())
    All_gradient_mat += gradient_matrix
    
    model.eval()
    with torch.no_grad():
        Y_train_pred = model(X_train)
        Y_test_pred = model(X_test)
        Y_train_plot.append(Y_train_pred.cpu().numpy())
        Y_test_plot.append(Y_test_pred.cpu().numpy())

        train_loss += criterion(Y_train_pred, Y_train).item()
        test_loss  += criterion(Y_test_pred, Y_test).item()

All_gradient_mat /= n_folds
print(f"Average Gradient Matrix:\n{All_gradient_mat}")


avg_train_loss = train_loss / n_folds
avg_test_loss = test_loss / n_folds
Loss_mat = [avg_train_loss, avg_test_loss]

print(f"Final Train Loss: {Loss_mat[0]:.4f}")
print(f"Final Test Loss:  {Loss_mat[1]:.4f}")



# SAVE FIG
Y_train_pred = np.concatenate(Y_train_plot, axis=0).reshape(-1, output_dim)
Y_test_pred  = np.concatenate(Y_test_plot, axis=0).reshape(-1, output_dim)
Y_train_true = np.concatenate(Y_train_true, axis=0).reshape(-1, output_dim)
Y_test_true  = np.concatenate(Y_test_true, axis=0).reshape(-1, output_dim)
# Combine all predictions and true values for plotting
all_vals = np.concatenate([
    Y_train_pred, Y_train_true,
    Y_test_pred,  Y_test_true
])



for output_idx in range(output_dim):
    plt.figure(figsize=(12, 6))
    min_val, max_val = all_vals.min(), all_vals.max()
    margin = 0.05 * (max_val - min_val)
    plot_min = min_val - margin
    plot_max = max_val + margin
    if test_D_flag == 1:
        title_train = f"Train: {output_name_list[output_idx]},All data"
        title_test = f"Test: {output_name_list[output_idx]},All data"
    else:
        title_train = f"Train: {output_name_list[output_idx]},Ggap:{D*735} nS"
        title_test = f"Test: {output_name_list[output_idx]},Ggap:{D*735} nS"
    plt.subplot(1, 2, 1)
    plt.scatter(Y_train_pred[:,output_idx], Y_train_true[:,output_idx], color='blue', alpha=0.5)
    plt.plot([Y_train_pred[:,output_idx].min(), Y_train_pred[:,output_idx].max()], [Y_train_pred[:,output_idx].min(), Y_train_pred[:,output_idx].max()], color='red')  

    plt.subplot(1, 2, 2)
    plt.scatter(Y_test_pred[:,output_idx], Y_test_true[:,output_idx], color='green', alpha=0.5)
    plt.plot([Y_test_pred[:,output_idx].min(), Y_test_pred[:,output_idx].max()], [Y_test_pred[:,output_idx].min(), Y_test_pred[:,output_idx].max()], color='red')  

    plt.tight_layout()
    plt.savefig(os.path.join(save_plot_path, f"KFold_result_{output_name_list[output_idx]}.png"))
    plt.close()

    print(f"{output_name_list[output_idx]} saved\n")

if test_D_flag != 1:
    grad_save_path = os.path.join(save_train_path, f"Grad_mat_D_{D}.txt")
    loss_save_path = os.path.join(save_train_path, f"Loss_D_{D}.txt")
else:
    grad_save_path = os.path.join(save_train_path, f"Grad_mat_D_All.txt")
    loss_save_path = os.path.join(save_train_path, f"Loss_D_All.txt")

np.savetxt(grad_save_path, All_gradient_mat.cpu().numpy(), fmt="%.5e")
np.savetxt(loss_save_path, Loss_mat, fmt="%.5e")

print(f"Gradient/ Loss matrix saved \n")

if test_D_flag == 1:
    visualize_feature_influence(All_gradient_mat, input_name_list, output_name_list, f"Feature Influence: {model_type} All data",save_plot_path)
else:
    visualize_feature_influence(All_gradient_mat, input_name_list, output_name_list, f"Feature Influence: {model_type} D={D}",save_plot_path)

end = time.time()
print(f"Total time taken: {end - start:.2f} seconds")