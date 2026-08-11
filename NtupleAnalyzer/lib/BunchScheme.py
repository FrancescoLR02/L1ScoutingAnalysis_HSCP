import requests
import urllib3

# Suppress insecure request warnings
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)

def get_cern_token(client_id, client_secret):
    """Fetches an OIDC token from CERN Keycloak using your App ID and Secret."""
    token_url = "https://auth.cern.ch/auth/realms/cern/protocol/openid-connect/token"
    payload = {
        "grant_type": "client_credentials",
        "client_id": client_id,
        "client_secret": client_secret,
        "audience": "cmsoms-prod"
    }
    
    response = requests.post(token_url, data=payload, verify=False)
    
    if response.status_code != 200:
        print(f"Error fetching token: {response.text}")
        response.raise_for_status()
        
    return response.json()["access_token"]

def generate_run_map(fill_numbers, token):
    """Fetches runs and bunch schemes using the Keycloak token."""
    base_url = "https://cmsoms.cern.ch/agg/api/v1"
    
    # Inject the token into the headers for all API requests
    headers = {
        "Authorization": f"Bearer {token}"
    }
    
    for fill in fill_numbers:
        print(f"   // {fill}")
        
        # 1. Fetch fill details
        fill_url = f"{base_url}/fills/{fill}"
        fill_resp = requests.get(fill_url, headers=headers, verify=False)
        
        if fill_resp.status_code != 200:
            print(f"   // Error fetching data for fill {fill}: Status {fill_resp.status_code}")
            continue
            
        fill_data = fill_resp.json()
        raw_bunch_scheme = fill_data['data']['attributes']['bunch_filling_scheme']
        
        # Apply your C++ alias naming logic
        if raw_bunch_scheme and "25ns" in raw_bunch_scheme:
             cpp_scheme_name = "scheme_25ns_v2"
        else:
             cpp_scheme_name = str(raw_bunch_scheme)
        
        # 2. Fetch runs for this fill
        runs_url = f"{base_url}/runs?filter[fill_number][EQ]={fill}&sort=-run_number"
        runs_resp = requests.get(runs_url, headers=headers, verify=False)
        runs_data = runs_resp.json()
        
        # 3. Print the output
        for run in runs_data.get('data', []):
            run_number = run['id']
            print(f'   {{{run_number}, "{cpp_scheme_name}"}},')

if __name__ == "__main__":
    # --- YOUR CREDENTIALS ---
    my_app_id = 'cms-l1scout-analysis'
    my_app_secret = 'XXX'  # Insert your actual secret here
    
    # --- YOUR FILLS ---
    fills_to_process = [11134, 11135] 
    
    try:
        # Step 1: Authenticate
        print("Authenticating with CERN Keycloak...")
        access_token = get_cern_token(my_app_id, my_app_secret)
        print("Authentication successful!\n")
        
        # Step 2: Generate Map
        generate_run_map(fills_to_process, access_token)
        
    except Exception as e:
        print(f"Script failed: {e}")