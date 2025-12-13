export type Json =
  | string
  | number
  | boolean
  | null
  | { [key: string]: Json | undefined }
  | Json[]

export type Database = {
  public: {
    Tables: {
      profiles: {
        Row: {
          id: string
          email: string | null
          full_name: string | null
          avatar_url: string | null
          created_at: string
          updated_at: string
        }
        Insert: {
          id: string
          email?: string | null
          full_name?: string | null
          avatar_url?: string | null
          created_at?: string
          updated_at?: string
        }
        Update: {
          id?: string
          email?: string | null
          full_name?: string | null
          avatar_url?: string | null
          created_at?: string
          updated_at?: string
        }
        Relationships: [
          {
            foreignKeyName: "profiles_id_fkey"
            columns: ["id"]
            isOneToOne: true
            referencedRelation: "users"
            referencedColumns: ["id"]
          }
        ]
      }
      user_api_keys: {
        Row: {
          id: string
          user_id: string
          provider: string
          encrypted_key: string
          created_at: string
          updated_at: string
        }
        Insert: {
          id?: string
          user_id: string
          provider: string
          encrypted_key: string
          created_at?: string
          updated_at?: string
        }
        Update: {
          id?: string
          user_id?: string
          provider?: string
          encrypted_key?: string
          created_at?: string
          updated_at?: string
        }
        Relationships: [
          {
            foreignKeyName: "user_api_keys_user_id_fkey"
            columns: ["user_id"]
            isOneToOne: false
            referencedRelation: "users"
            referencedColumns: ["id"]
          }
        ]
      }
      analysis_results: {
        Row: {
          id: string
          user_id: string
          analysis_type: string
          title: string | null
          description: string | null
          input_data: Json | null
          results: Json | null
          settings: Json | null
          created_at: string
          updated_at: string
        }
        Insert: {
          id?: string
          user_id: string
          analysis_type: string
          title?: string | null
          description?: string | null
          input_data?: Json | null
          results?: Json | null
          settings?: Json | null
          created_at?: string
          updated_at?: string
        }
        Update: {
          id?: string
          user_id?: string
          analysis_type?: string
          title?: string | null
          description?: string | null
          input_data?: Json | null
          results?: Json | null
          settings?: Json | null
          created_at?: string
          updated_at?: string
        }
        Relationships: [
          {
            foreignKeyName: "analysis_results_user_id_fkey"
            columns: ["user_id"]
            isOneToOne: false
            referencedRelation: "users"
            referencedColumns: ["id"]
          }
        ]
      }
      analysis_sessions: {
        Row: {
          id: string
          user_id: string
          session_name: string | null
          analysis_type: string | null
          status: string
          progress: Json | null
          created_at: string
          updated_at: string
        }
        Insert: {
          id?: string
          user_id: string
          session_name?: string | null
          analysis_type?: string | null
          status?: string
          progress?: Json | null
          created_at?: string
          updated_at?: string
        }
        Update: {
          id?: string
          user_id?: string
          session_name?: string | null
          analysis_type?: string | null
          status?: string
          progress?: Json | null
          created_at?: string
          updated_at?: string
        }
        Relationships: [
          {
            foreignKeyName: "analysis_sessions_user_id_fkey"
            columns: ["user_id"]
            isOneToOne: false
            referencedRelation: "users"
            referencedColumns: ["id"]
          }
        ]
      }
      feedback_comments: {
        Row: {
          id: string
          user_id: string
          comment_type: 'feature_request' | 'bug_report'
          name: string | null
          email: string | null
          message: string
          is_approved: boolean
          approval_token: string
          created_at: string
          updated_at: string
        }
        Insert: {
          id?: string
          user_id: string
          comment_type: 'feature_request' | 'bug_report'
          name?: string | null
          email?: string | null
          message: string
          is_approved?: boolean
          approval_token?: string
          created_at?: string
          updated_at?: string
        }
        Update: {
          id?: string
          user_id?: string
          comment_type?: 'feature_request' | 'bug_report'
          name?: string | null
          email?: string | null
          message?: string
          is_approved?: boolean
          approval_token?: string
          created_at?: string
          updated_at?: string
        }
        Relationships: [
          {
            foreignKeyName: "feedback_comments_user_id_fkey"
            columns: ["user_id"]
            isOneToOne: false
            referencedRelation: "users"
            referencedColumns: ["id"]
          }
        ]
      }
    }
    Views: {
      [_ in never]: never
    }
    Functions: {
      [_ in never]: never
    }
    Enums: {
      [_ in never]: never
    }
    CompositeTypes: {
      [_ in never]: never
    }
  }
}

export type Provider = 'openrouter' | 'anthropic' | 'openai'

export type AnalysisType = 'batch' | 'symphony' | 'scoring' | 'subclustering' | 'annotation-boost'

export type SessionStatus = 'active' | 'completed' | 'archived'

export type Profile = Database['public']['Tables']['profiles']['Row']
export type UserApiKey = Database['public']['Tables']['user_api_keys']['Row']
export type AnalysisResult = Database['public']['Tables']['analysis_results']['Row']
export type AnalysisSession = Database['public']['Tables']['analysis_sessions']['Row']
export type FeedbackComment = Database['public']['Tables']['feedback_comments']['Row']
export type CommentType = 'feature_request' | 'bug_report'