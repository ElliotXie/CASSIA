import os
from flask import Flask, render_template, request, jsonify, Response, session
import pandas as pd
from werkzeug.utils import secure_filename
from version2_nohuman_wrapper_926 import run_cell_type_analysis
from openai import OpenAI
import json
import io
import webbrowser
from threading import Timer

# Import the formatting function
from processcode import format_annotation_output

def create_app():
    app = Flask(__name__)
    app.secret_key = 'your_secret_key_here'  # Replace with a real secret key

    ALLOWED_EXTENSIONS = {'csv'}

    def allowed_file(filename):
        return '.' in filename and filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS

    @app.route('/', methods=['GET', 'POST'])
    def index():
        if request.method == 'POST':
            # Handle form submission
            openai_key = request.form['openai_key']
            tissue = request.form['tissue']
            species = request.form['species']
            
            # Handle file upload
            if 'file' not in request.files:
                return jsonify({'error': 'No file part'})
            file = request.files['file']
            if file.filename == '':
                return jsonify({'error': 'No selected file'})
            if file and allowed_file(file.filename):
                # Read the CSV file directly from the request
                df = pd.read_csv(file)
                
                # Store the DataFrame in the session
                session['df'] = df.to_json()
                session['openai_key'] = openai_key
                session['tissue'] = tissue
                session['species'] = species
                
                return render_template('results.html', total_rows=len(df))
            else:
                return jsonify({'error': 'Invalid file type'})
        
        return render_template('index.html')

    @app.route('/stream_results')
    def stream_results():
        openai_key = session.get('openai_key')
        tissue = session.get('tissue')
        species = session.get('species')
        df = pd.read_json(session.get('df'))
        
        client = OpenAI(api_key=openai_key)

        def generate():
            yield f"data: {json.dumps({'type': 'console', 'message': 'Starting analysis...'})}\n\n"
            
            for index, row in df.head(2).iterrows():
                broad_cell_type = row['Broad cell type']
                marker_list = row['Top 10 Markers'].split(', ')
                
                yield f"data: {json.dumps({'type': 'console', 'message': f'Starting final annotation for {broad_cell_type}...'})}\n\n"
                
                result, conversation_history = run_cell_type_analysis("gpt-4", 0, marker_list, tissue, species, "no")
                
                yield f"data: {json.dumps({'type': 'console', 'message': f'Validating annotation for {broad_cell_type}...'})}\n\n"
                
                if result:
                    yield f"data: {json.dumps({'type': 'console', 'message': f'Formatting final results for {broad_cell_type}...'})}\n\n"
                    
                    # Format the conversation history
                    formatted_history = []
                    for role, message in conversation_history:
                        try:
                            if role == "Final Annotation Agent":
                                formatted_message = format_annotation_output(message)
                            else:
                                formatted_message = message
                            formatted_history.append((role, formatted_message))
                        except Exception as e:
                            print(f"Error formatting message: {e}")
                            formatted_history.append((role, message))
                    
                    data = {
                        "type": "result",
                        "broad_cell_type": broad_cell_type,
                        "analysis_result": result,
                        "conversation_history": formatted_history,
                        "confidence_score": result.get("confidence_score")
                    }
                    yield f"data: {json.dumps(data)}\n\n"

        return Response(generate(), content_type='text/event-stream')

    @app.route('/about')
    def about():
        return render_template('about.html')

    return app

def run_cell_type_analysis_app(port=5000):
    app = create_app()
    
    # Open the browser after a short delay
    def open_browser():
        webbrowser.open(f'http://127.0.0.1:{port}/')
    
    Timer(1.5, open_browser).start()
    
    # Run the Flask app
    app.run(port=port, debug=True, use_reloader=False)

if __name__ == '__main__':
    run_cell_type_analysis_app()